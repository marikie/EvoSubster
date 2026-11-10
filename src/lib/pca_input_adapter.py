"""Normalize local legacy PCA inputs without weakening v2 manifest validation."""
import argparse
import json
import os
from pathlib import Path
import re

import run_manifest as runs
from train_identity import IDENTITY_SOURCE, parse_identity

TRIO = re.compile(r'[^_]+1_[^_]+2_[^_]+3')
DATE = re.compile(r'(?:[0-9]{8}|[0-9]{14})')


def raw_manifest(run):
    path = run / runs.MANIFEST
    value = json.loads(path.read_text())
    if not isinstance(value, dict):
        raise ValueError(f'Manifest must be an object: {path}')
    if 'schema_version' in value and value['schema_version'] != runs.SCHEMA:
        raise ValueError(f'Unsupported manifest schema: {path}')
    return value


def ignored(name):
    analysis_dirs = ('spectrum_shape_bench_noncoding', 'spectrum_shape_bench_whole')
    return name in ('archive', 'tmp', 'test', *analysis_dirs) or name.startswith(('pca', 'bench_', '_'))


def discover(root):
    root = Path(root).resolve()
    if (root / runs.MANIFEST).is_file():
        value = raw_manifest(root)
        if 'schema_version' in value:
            value = runs.load(root)
            return [str(root)] if value['status'] == 'complete' else []
        return [str(root)]
    legacy = []
    modern = {}
    visited = set()
    for current, directories, _ in os.walk(root, followlinks=True):
        current = Path(current)
        actual = current.resolve()
        if actual in visited:
            directories[:] = []
            continue
        visited.add(actual)
        directories[:] = sorted(d for d in directories if not ignored(d))
        if TRIO.fullmatch(current.name) and (current == root or current.parent == root):
            dated = sorted(p for p in current.iterdir() if p.is_dir() and DATE.fullmatch(p.name))
            if dated:
                latest = dated[-1]
                value = raw_manifest(latest)
                if 'schema_version' not in value:
                    legacy.append(str(latest.resolve()))
                    directories[:] = []
                    continue
        if (current / runs.MANIFEST).is_file():
            value = raw_manifest(current)
            if 'schema_version' in value:
                value = runs.load(current)
                if value['status'] == 'complete':
                    key = value['trio_key']
                    if key not in modern or value['created_at'] > modern[key][0]:
                        modern[key] = value['created_at'], str(actual)
                directories[:] = []
    result = legacy + [run for _, run in sorted(modern.values())]
    if not result:
        raise ValueError(f'No supported completed PCA runs found: {root}')
    return sorted(set(result))


def safe_component(value, field):
    if not isinstance(value, str) or not value or Path(value).name != value or value in ('.', '..'):
        raise ValueError(f'Invalid legacy {field}: {value}')
    return value


def read(run):
    run = Path(run).resolve()
    value = raw_manifest(run)
    if 'schema_version' in value:
        return runs.validated(run, complete=True)
    if not DATE.fullmatch(run.name) or not TRIO.fullmatch(run.parent.name):
        raise ValueError(f'Legacy PCA requires a dated run in a trio directory: {run}')
    organisms = value.get('organisms')
    if not isinstance(organisms, list) or len(organisms) != 3 or not all(isinstance(e, dict) for e in organisms):
        raise ValueError(f'Invalid legacy organisms: {run}')
    slots = {entry.get('slot'): dict(entry) for entry in organisms}
    if set(slots) != {'org1', 'org2', 'org3'}:
        raise ValueError(f'Invalid legacy slots: {run}')
    artifacts = []
    issues = []
    for slot, entry in slots.items():
        short = safe_component(entry.get('short_name'), 'short_name')
        accession = safe_component(entry.get('accession'), 'accession')
        entry['role'] = 'outgroup' if slot == 'org1' else 'ingroup'
        entry['metadata_json'] = f'metadata/{short}_{accession}.json'
        taxonomy = f'metadata/taxonomy_{short}_{accession}.json'
        entry['taxonomy_json'] = taxonomy if (run / taxonomy).is_file() else None
        if slot == 'org1':
            continue
        for kind, suffix in (('single', ''), ('single_ncds', '_ncds')):
            filename = f'{accession}_{short}_{run.name}{suffix}.tsv'
            matches = [p for p in (run / filename, run / 'statistics' / short / 'singlenuc' / filename) if p.is_file()]
            if len(matches) > 1:
                raise ValueError(f'Ambiguous canonical TSV for {slot}: {filename}')
            if matches:
                artifacts.append(dict(kind=kind, slot=slot, path=str(matches[0].relative_to(run))))
    for ref, query in ((1, 2), (1, 3), (2, 3)):
        filename = f"{slots[f'org{ref}']['short_name']}2{slots[f'org{query}']['short_name']}_{run.name}.train"
        matches = [p for p in (run / filename, run / 'intermediateFiles' / filename) if p.is_file()]
        if len(matches) > 1:
            raise ValueError(f'Ambiguous LAST training files: {filename}')
        if not matches:
            issues.append(f'Missing LAST identity: {filename}')
            continue
        try:
            identity = parse_identity(matches[0])
        except (ValueError, OSError) as error:
            issues.append(str(error))
            continue
        artifacts.append(dict(kind='train', slot=f'org{ref}-org{query}',
                              identity_percent=identity, identity_source=IDENTITY_SOURCE, path=str(matches[0].relative_to(run))))
    return dict(identity_source=IDENTITY_SOURCE, layout='legacy', trio_key=run.parent.name, date=run.name,
                organisms=list(slots.values()), artifacts=artifacts, input_issues=issues)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['discover', 'read'])
    parser.add_argument('path')
    args = parser.parse_args()
    print(json.dumps(discover(args.path) if args.command == 'discover' else read(args.path)))
