#!/usr/bin/env python3
"""Create a checked collection of modern runs and legacy trio directories."""
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
import run_manifest as runs
from pca_identity_filter import evaluate
from train_identity import IDENTITY_SOURCE, parse_identity

TRIO = re.compile(r'[^_]+1_[^_]+2_[^_]+3')
DATE = re.compile(r'(?:[0-9]{8}|[0-9]{14})')
OPTIONS = ('idt_threshold', 'ingroup_idt_threshold', 'outgroup_idt_threshold',
           'min_idt_gap_between_out_in')
COLUMNS = ('trio_id', 'layout', 'run_path', 'link_source', 'idt_12', 'idt_13',
           'idt_23', 'gap_12', 'gap_13', 'genus_org1', 'genus_org2', 'genus_org3',
           'genus_condition', 'gap_waived', 'ingroup_idt_threshold',
           'outgroup_idt_threshold', 'min_idt_gap_between_out_in', 'included', 'reason',
           'link_action', 'identity_source')


def read_manifest(run):
    path = run / runs.MANIFEST
    try:
        manifest = json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise ValueError(f'Cannot read manifest {path}: {error}') from error
    if not isinstance(manifest, dict):
        raise ValueError(f'Manifest must be an object: {path}')
    if 'schema_version' in manifest and manifest['schema_version'] != runs.SCHEMA:
        raise ValueError(f'Unsupported manifest schema: {path}')
    return manifest


def legacy_candidate(trio, output):
    dated = sorted((p for p in trio.iterdir()
                    if p.is_dir() and DATE.fullmatch(p.name)
                    and not p.resolve().is_relative_to(output)), key=lambda p: p.name)
    if not dated:
        raise ValueError(f'Legacy trio has no usable dated run: {trio}')
    run = dated[-1].resolve()
    manifest = read_manifest(run)
    if 'schema_version' in manifest:
        return None
    organisms = manifest.get('organisms')
    if not isinstance(organisms, list) or sorted(e.get('slot', '') for e in organisms if isinstance(e, dict)) != ['org1', 'org2', 'org3'] or len(organisms) != 3:
        raise ValueError(f'Invalid legacy organism slots: {run}')
    slots = {entry['slot']: entry for entry in organisms}
    names = [slots[slot].get('short_name') for slot in ('org1', 'org2', 'org3')]
    if any(not isinstance(name, str) or not name or Path(name).name != name or name in ('.', '..') for name in names):
        raise ValueError(f'Invalid legacy short names: {run}')
    identities = {}
    issues = []
    for ref, query in ((1, 2), (1, 3), (2, 3)):
        filename = f'{names[ref - 1]}2{names[query - 1]}_{run.name}.train'
        files = [p for p in (run / filename, run / 'intermediateFiles' / filename) if p.is_file()]
        if len(files) > 1:
            raise ValueError(f'Ambiguous training files for {trio.name}: {filename}')
        key = f'idt_{ref}{query}'
        if not files:
            issues.append(f'{key}: missing training file {filename}')
            continue
        try:
            identities[key] = parse_identity(files[0])
        except (ValueError, OSError) as error:
            issues.append(f'{key}: {error}')
    return dict(key=trio.name, source=trio.resolve(), run=run, layout='legacy',
                organisms=slots, identities=identities, issues=issues)


def modern_candidate(run, manifest):
    identities = {}
    for artifact in manifest['artifacts']:
        if artifact['kind'] == 'train':
            pair = artifact['slot'].replace('org', '').replace('-', '')
            identities['idt_' + pair] = artifact['identity_percent']
    return dict(
        key=manifest['trio_key'], source=run, run=run, layout='v2',
        organisms={entry['slot']: entry for entry in manifest['organisms']},
        identities=identities, issues=[],
    )


def discover(root, output):
    """Legacy discovery is shallow; modern manifests retain recursive discovery."""
    if (root / runs.MANIFEST).is_file():
        manifest = read_manifest(root)
        if 'schema_version' in manifest:
            manifest = runs.load(root)
            if manifest['status'] != 'complete':
                return []
            manifest = runs.validated(root, complete=True)
            return [modern_candidate(root, manifest)]
    trios = [root] if TRIO.fullmatch(root.name) else [
        p for p in sorted(root.iterdir())
        if p.is_dir() and TRIO.fullmatch(p.name) and not p.resolve().is_relative_to(output)
    ]
    modern = []
    modern_runs = set()
    visited = set()
    for current, directories, files in os.walk(root, followlinks=True):
        actual = Path(current).resolve()
        if actual.is_relative_to(output) or actual in visited:
            directories[:] = []
            continue
        visited.add(actual)
        directories[:] = sorted(d for d in directories
                                if not (Path(current) / d).resolve().is_relative_to(output))
        if Path(current).name == 'metadata' and 'metadata_manifest.json' in files:
            run = Path(current).parent.resolve()
            manifest = read_manifest(run)
            if 'schema_version' not in manifest:
                continue
            manifest = runs.load(run)
            modern_runs.add(run)
            if manifest['status'] == 'complete':
                modern.append((manifest['trio_key'], manifest['created_at'], run))
    latest = {}
    for key, timestamp, run in sorted(modern):
        if key not in latest or timestamp > latest[key][0]:
            latest[key] = timestamp, run
    legacy = []
    for trio in trios:
        # Modern run roots with trio-like names do not require dated children.
        if (trio / runs.MANIFEST).is_file():
            manifest = read_manifest(trio)
            if 'schema_version' in manifest:
                continue
        has_dates = any(p.is_dir() and DATE.fullmatch(p.name) for p in trio.iterdir())
        if not has_dates and any(run.is_relative_to(trio.resolve()) for run in modern_runs):
            continue
        candidate = legacy_candidate(trio, output)
        if candidate is not None:
            legacy.append(candidate)
    candidates = legacy
    for key, (_, run) in sorted(latest.items()):
        manifest = runs.validated(run, complete=True)
        candidates.append(modern_candidate(run, manifest))
    if not candidates and not modern_runs:
        raise ValueError(f'No completed v2 runs or legacy trio directories found: {root}')
    return candidates


def plan_links(roots, output, options):
    if any(output == root for root in roots):
        raise ValueError('Output cannot also be an input')
    if output.is_symlink() or (output.exists() and not output.is_dir()):
        raise ValueError(f'Output must be a directory: {output}')
    planned = {}
    for root in roots:
        if not root.is_dir():
            raise ValueError(f'Input is not a directory: {root}')
        for candidate in discover(root, output):
            key = candidate['key']
            if not isinstance(key, str) or not key or Path(key).name != key or key in ('.', '..', 'link_run_settings.json', 'link_trio_filter.csv'):
                raise ValueError(f'Invalid trio key: {key}')
            if key in planned:
                previous = planned[key]
                if previous['source'] != candidate['source'] or previous['run'] != candidate['run']:
                    raise ValueError(f'Conflicting runs for trio: {key}')
                continue
            decision = evaluate(dict(options=options, identities=candidate['identities'],
                                     organisms=candidate['organisms']))
            genera = decision.pop('genera')
            decision.pop('thresholds')
            decision.update({f'genus_{slot}': genera.get(slot) for slot in ('org1', 'org2', 'org3')})
            decision.update(trio_id=key, layout=candidate['layout'],
                            run_path=str(candidate['run']), link_source=str(candidate['source']))
            if candidate['issues']:
                decision['reason'] = '; '.join(filter(None, [decision['reason'], *candidate['issues']]))
            candidate['decision'] = decision
            planned[key] = candidate
    for key, candidate in planned.items():
        target = output / key
        if target.exists() or target.is_symlink():
            if not target.is_symlink() or target.resolve() != candidate['source']:
                raise ValueError(f'Destination collision: {target}')
        if candidate['decision']['included']:
            action = 'reused' if target.is_symlink() else 'created'
        else:
            action = 'kept_existing' if target.is_symlink() else 'skipped'
        candidate['decision']['link_action'] = action
    for name in ('link_run_settings.json', 'link_trio_filter.csv'):
        target = output / name
        if target.is_symlink() or (target.exists() and not target.is_file()):
            raise ValueError(f'Audit destination collision: {target}')
    return planned


def main():
    parser = argparse.ArgumentParser(description='Link completed v2 runs and legacy trios with optional identity thresholds')
    parser.add_argument('--output-dir', required=True)
    descriptions = {
        'idt_threshold': 'Common strict identity lower bound; individual bounds override it',
        'ingroup_idt_threshold': 'Strict ingroup identity lower bound',
        'outgroup_idt_threshold': 'Strict lower bound for both outgroup-ingroup identities',
        'min_idt_gap_between_out_in': 'Strict gap in percentage points; waived for same-genus ingroups with a different-genus outgroup',
    }
    for option, description in descriptions.items():
        parser.add_argument('--' + option.replace('_', '-'), type=float,
                            help=description + ' (finite [0,100]; disabled unless specified)')
    parser.add_argument('input_dirs', nargs='+')
    args = parser.parse_args()
    options = {key: getattr(args, key) for key in OPTIONS}
    effective = evaluate(dict(options=options))['thresholds']
    roots = list(dict.fromkeys(Path(value).resolve() for value in args.input_dirs))
    # Check symlink output before resolving so audit writes cannot redirect.
    destination = Path(args.output_dir).absolute()
    if destination.is_symlink():
        raise ValueError(f'Output must not be a symlink: {destination}')
    output = destination.resolve()
    planned = plan_links(roots, output, options)
    included = {key for key, candidate in planned.items() if candidate['decision']['included']}
    preserved_existing = sum(
        target.is_symlink() and target.name not in included
        for target in output.iterdir()
    ) if output.is_dir() else 0
    created = reused = excluded = 0
    output.mkdir(parents=True, exist_ok=True)
    for key, candidate in sorted(planned.items()):
        if not candidate['decision']['included']:
            excluded += 1
            continue
        target = output / key
        if target.is_symlink():
            reused += 1
        else:
            target.symlink_to(os.path.relpath(candidate['source'], output), target_is_directory=True)
            created += 1
    settings = dict(
        arguments=sys.argv[1:], requested_thresholds=options, effective_thresholds=effective,
        disabled_conditions=[key for key, value in effective.items() if value is None],
        identity_source=IDENTITY_SOURCE,
        comparison='>', identity_unit='percent', gap_unit='percentage points',
        genus_waiver='same-genus ingroups and different-genus outgroup waive only gap',
        input_dirs=[str(root) for root in roots], output_dir=str(output),
        recorded_at=datetime.now(timezone.utc).isoformat(), trios=len(planned),
        created=created, reused=reused, excluded=excluded, included=created + reused,
        preserved_existing=preserved_existing,
    )
    runs.atomic_json(output / 'link_run_settings.json', settings)
    with (output / 'link_trio_filter.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(candidate['decision'] for _, candidate in sorted(planned.items()))
    print(f'Created: {created}; reused: {reused}; excluded: {excluded}; '
          f'preserved_existing: {preserved_existing}')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as error:
        sys.exit(f'Error: {error}')
