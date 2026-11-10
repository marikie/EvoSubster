#!/usr/bin/env python3
"""Explicit, portable completion records for independently reusable pipeline stages.

Only register after the producing command has succeeded. Existing files alone never
establish completion. Source runs are read under shared locks and never modified.
"""
import argparse
import copy
from contextlib import contextmanager, ExitStack
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

import run_manifest as manifest
from train_identity import parse_identity

ROOT = Path(__file__).resolve().parents[2]
SLOTS = ('org1', 'org2', 'org3')
TRAINS = ('org1-org2', 'org1-org3', 'org2-org3')
PAIRS = TRAINS[:2]
POLICIES = ('same', 'allow-existing-pairs')
ERRORS = (ValueError, OSError, KeyError, TypeError, AttributeError)


def _check(kind, slot):
    allowed = SLOTS if kind in ('gc', 'coding') else TRAINS if kind == 'train' else PAIRS if kind == 'pair' else ()
    if slot not in allowed:
        raise ValueError(f'Unsupported reusable stage: {kind}:{slot}')


def _implementations(kind):
    align = Path(os.environ.get('ALIGN_DIR_OVERRIDE', ROOT/'src/align'))
    metrics = Path(os.environ.get('METRICS_DIR_OVERRIDE', ROOT/'src/metrics'))
    paths = {'gc': [metrics/'gc_content.sh'],
             'coding': [metrics/'coding_region_percentage.py'],
             'train': [align/'last_train.sh', ROOT/'src/lib/last_options.py'],
             'pair': [align/'one2one.sh', align/'last_train.sh', ROOT/'src/lib/last_options.py']}[kind]
    return {p.name: manifest.sha256(p) for p in paths}


def _without_threads(words):
    result, index = [], 0
    while index < len(words):
        word = words[index]
        if word == '-P':
            if index + 1 >= len(words):
                raise ValueError('Missing LAST thread count')
            index += 2
        elif re.fullmatch(r'-P[0-9]+', word):
            index += 1
        else:
            result.append(word)
            index += 1
    return result


def _signature(value, kind, slot, train_sha=None):
    _check(kind, slot)
    organisms = {o['slot']: o for o in value['organisms']}
    options = value['options']
    inputs, settings = {}, {}
    if kind in ('gc', 'coding'):
        inputs['fasta'] = organisms[slot]['fasta_sha256']
        if kind == 'coding':
            inputs['gff'] = organisms[slot]['gff_sha256']
            if inputs['gff'] is None:
                raise ValueError('Coding stage requires a GFF')
    else:
        ref, query = slot.split('-')
        inputs = {'reference': organisms[ref]['fasta_sha256'], 'query': organisms[query]['fasta_sha256']}
        last = options['last']
        if last.get('schema') != 1:
            raise ValueError('Unsupported LAST option schema')
        settings = {k: copy.deepcopy(last[k]) for k in ('lastdb_train', 'last_train')}
        if kind == 'pair':
            inputs['train'] = train_sha
            settings.update({k: copy.deepcopy(last[k]) for k in ('lastdb_align', 'last_split')})
            settings['lastal'] = _without_threads(last['lastal'])
        from last_options import scientific_option_files
        external=scientific_option_files(last,keys=settings.keys())
        if external:inputs['option_files']=external
    return dict(schema=1, kind=kind, inputs=inputs, options=settings,
                algorithm=options.get('algorithm', value.get('algorithm', manifest.ALGORITHM)),
                implementation=_implementations(kind))


def _path(run, relative):
    if not isinstance(relative, str) or Path(relative).is_absolute() or '..' in Path(relative).parts:
        raise ValueError(f'Output must be a relative path inside run: {relative}')
    run = Path(run).resolve()
    raw = run/relative
    # Reject every symlink component, including links that remain inside the run.
    current = run
    for part in Path(relative).parts:
        current = current/part
        if current.is_symlink():
            raise ValueError(f'Symlink output is not reusable: {relative}')
    return manifest.artifact_path(run, relative)


def validate_maf(path):
    """Validate two sequence rows per block, aligned lengths and coordinates."""
    blocks, rows, active = 0, [], False
    def finish():
        if active and (len(rows) != 2 or len(rows[0]) != len(rows[1])):
            raise ValueError(f'Incomplete pairwise MAF block: {path}')
    with Path(path).open() as stream:
        for number, line in enumerate(stream, 1):
            fields = line.split()
            if not fields or line.lstrip().startswith('#'):
                continue
            if fields[0] == 'a':
                finish()
                blocks += 1
                rows, active = [], True
            elif fields[0] == 's':
                if not active or len(fields) != 7:
                    raise ValueError(f'Invalid MAF sequence row: {path}:{number}')
                _, name, start, size, strand, total, sequence = fields
                start, size, total = int(start), int(size), int(total)
                if (start < 0 or size <= 0 or total <= 0 or start+size > total
                        or strand not in ('+', '-') or len(sequence.replace('-', '')) != size
                        or not re.fullmatch('[A-Za-z*-]+', sequence)):
                    raise ValueError(f'Invalid MAF coordinate or length: {path}:{number}')
                rows.append(sequence)
                if len(rows) > 2:
                    raise ValueError(f'MAF is not pairwise: {path}')
            elif fields[0] not in ('q', 'p', 'i') or not active:
                raise ValueError(f'Unsupported MAF row: {path}:{number}')
    finish()
    if not blocks:
        raise ValueError(f'Empty pairwise MAF: {path}')
    return blocks


def _outputs(run, kind, outputs):
    outputs = dict(outputs)
    if kind in ('gc', 'coding'):
        if len(outputs) != 1 or not set(outputs).issubset({'result', kind}):
            raise ValueError('Metric stage requires result=PATH')
        outputs = {'result': next(iter(outputs.values()))}
    elif kind == 'train' and set(outputs) != {'train'}:
        raise ValueError('Training stage requires train=PATH')
    elif kind == 'pair' and set(outputs) not in ({'train', 'maf'}, {'train', 'maf_sorted'}):
        raise ValueError('Pair stage requires train=PATH and maf=PATH (or maf_sorted=PATH)')
    result = []
    for role, relative in outputs.items():
        path = _path(run, relative)
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f'Missing or empty output: {path}')
        if role == 'train':
            try:parse_identity(path)
            except ValueError as error:
                raise ValueError(f'incomplete training output: {path}: {error}') from error
        elif role.startswith('maf'):
            validate_maf(path)
        result.append(dict(role=role, path=relative, sha256=manifest.sha256(path)))
    return result


def _file(record, role):
    return next(o for o in record['outputs'] if o['role'] == role)


def _verify_outputs(run, record):
    roles = {o['role']: o['path'] for o in record['outputs']}
    if len(roles) != len(record['outputs']):
        raise ValueError('Duplicate output roles')
    checked = _outputs(run, record['kind'], roles)
    if checked != record['outputs']:
        raise ValueError('Stage output checksum mismatch')


def _compatible(value, kind, slot, record, allow_old=False):
    if record.get('schema') != 1 or record.get('kind') != kind:
        return False
    sha = _file(record, 'train')['sha256'] if kind == 'pair' else None
    if record['signature'] != _signature(value, kind, slot, sha):
        return False
    if kind in ('train', 'pair') and not allow_old:
        return record['producer']['last_version'] == value['options']['last_version']
    return True


def _verified(run, value, kind, slot, allow_old=False):
    record = value.get('stage_records', {}).get(f'{kind}:{slot}')
    if not record or record.get('slot') != slot or not _compatible(value, kind, slot, record, allow_old):
        raise ValueError('Missing completion record or stage signature mismatch')
    _verify_outputs(run, record)
    if kind == 'pair':
        training = value.get('stage_records', {}).get(f'train:{slot}')
        if not training or not _compatible(value, 'train', slot, training, allow_old):
            raise ValueError('Missing matching training completion record')
        _verify_outputs(run, training)
        if _file(training, 'train')['sha256'] != _file(record, 'train')['sha256']:
            raise ValueError('Pair training profile disagrees with completed training')
        if training['producer']['last_version'] != record['producer']['last_version']:
            raise ValueError('Pair and training producer versions disagree')
    return record


def _ready(run, value, kind, slot):
    try:
        record = value.get('stage_records', {}).get(f'{kind}:{slot}', {})
        old = kind == 'pair' and record.get('provenance', {}).get('policy') == 'allow-existing-pairs'
        _verified(run, value, kind, slot, old)
        return True
    except ERRORS:
        if kind == 'train' and slot in PAIRS:
            try:
                pair = value.get('stage_records', {}).get(f'pair:{slot}', {})
                if pair.get('provenance', {}).get('policy') != 'allow-existing-pairs':
                    return False
                _verified(run, value, 'pair', slot, True)
                return True
            except ERRORS:
                pass
        return False


def ready(run, kind, slot):
    _check(kind, slot)
    try:
        _path(run, manifest.MANIFEST)
        return _ready(run, manifest.load(run), kind, slot)
    except ERRORS:
        return False


def reconcile(run, value):
    """Invalidate stale native stages before existing files can bypass execution.

    Called by prepare while holding the run lock. Returning fewer artifact entries
    lets prepare remove only the generated files whose completion is no longer valid.
    """
    records=value.get('stage_records',{})
    invalid={key for key,r in records.items() if not _ready(run,value,r['kind'],r['slot'])}
    if not invalid:
        return value
    discarded={o['path'] for key in invalid for o in records[key]['outputs']}
    intact={o['path'] for key,r in records.items() if key not in invalid for o in r['outputs']}
    discarded-=intact
    downstream=any(records[key]['kind'] in ('train','pair') for key in invalid)
    kept=[]
    for artifact in value['artifacts']:
        name=artifact['path']
        if name in discarded or any(name==p+suffix for p in discarded for suffix in ('.complete.json','.last_options.json')):
            continue
        if downstream:
            if artifact['kind'] not in ('gc','coding','train','auxiliary'):
                continue
            if artifact['kind']=='auxiliary' and name.startswith(('intermediateFiles/','figs/')):
                if name not in intact and not any(name==p+'.complete.json' for p in intact):
                    continue
        kept.append(artifact)
    value['artifacts']=kept
    value['stage_records']={key:r for key,r in records.items() if key not in invalid}
    value['status']='pending'
    return value


def _artifact(value, run, kind, slot, output):
    entry = dict(kind=kind, slot=slot, path=output['path'], sha256=output['sha256'])
    if kind == 'train':
        entry['identity_percent'] = parse_identity(_path(run, output['path']))
    # Every path and every logical artifact appears at most once.
    value['artifacts'] = [a for a in value['artifacts'] if a['path'] != entry['path']
                          and (kind == 'auxiliary' or (a['kind'], a['slot']) != (kind, slot))]
    value['artifacts'].append(entry)


def _record_artifacts(value, run, record):
    for output in record['outputs']:
        role = output['role']
        kind = 'train' if role == 'train' else record['kind'] if role == 'result' else 'auxiliary'
        _artifact(value, run, kind, record['slot'] if kind != 'auxiliary' else None, output)
    # Preserve existing, verified sidecars so run preparation does not remove them.
    for output in record['outputs']:
        path = _path(run, output['path'])
        candidates = [Path(str(path)+'.complete.json')]
        if output['role'].startswith('maf'):
            candidates.append(Path(str(path)+'.last_options.json'))
        for sidecar in candidates:
            if not sidecar.is_file():
                continue
            data = json.loads(sidecar.read_text())
            if sidecar.name.endswith('.last_options.json'):
                from last_options import alignment_settings
                valid = (data == value['options']['last']
                         or data == alignment_settings(value['options']['last']))
            else:
                valid = data.get('sha256') == output['sha256']
                if output['role'] == 'train':
                    valid = data == {'sha256': output['sha256']}
            if valid:
                relative = sidecar.relative_to(Path(run).resolve()).as_posix()
                _artifact(value, run, 'auxiliary', None, dict(path=relative, sha256=manifest.sha256(sidecar)))


def register(run, kind, slot, outputs):
    """Register a caller-confirmed successful native stage; never infer completion."""
    _check(kind, slot)
    run = Path(run).resolve()
    _path(run, manifest.MANIFEST)
    value = manifest.load(run)
    checked = _outputs(run, kind, outputs)
    old = value.get('stage_records', {}).get(f'{kind}:{slot}')
    if old and _ready(run, value, kind, slot) and old['outputs'] == checked:
        return old
    train_sha = next((o['sha256'] for o in checked if o['role'] == 'train'), None)
    record = dict(schema=1, kind=kind, slot=slot, signature=_signature(value, kind, slot, train_sha),
                  producer=dict(last_version=value['options'].get('last_version'), implementation=_implementations(kind)),
                  outputs=checked)
    value.setdefault('stage_records', {})[f'{kind}:{slot}'] = record
    if kind == 'pair':
        _verified(run, value, kind, slot)
    _record_artifacts(value, run, record)
    manifest.atomic_json(run/manifest.MANIFEST, value)
    return record


def _discover(source, target=None):
    root = Path(source).resolve()
    if not root.is_dir():
        raise ValueError(f'Reuse source is not a directory: {root}')
    explicit = (root/manifest.MANIFEST).is_file()
    candidates = [root] if explicit else []
    if not explicit:
        count = 0
        for current, directories, files in os.walk(root, followlinks=False):
            count += 1
            if count > 100000:
                raise ValueError('Reuse source discovery exceeded directory limit')
            here = Path(current)
            directories[:] = [d for d in directories if not (here/d).is_symlink()]
            if len(here.relative_to(root).parts) >= 12:
                directories[:] = []
            if here.name == 'metadata' and 'metadata_manifest.json' in files:
                candidates.append(here.parent.resolve())
    result = []
    for run in sorted(set(candidates)):
        if not run.is_relative_to(root):
            raise ValueError('Discovered run escapes reuse source root')
        if target is not None:
            target = Path(target).resolve()
            overlap = target == run or target.is_relative_to(run) or run.is_relative_to(target)
            if overlap:
                if explicit:
                    raise ValueError('Reuse source and target must not overlap')
                continue
        result.append(run)
    return result, explicit


@contextmanager
def _sources(source, target=None):
    candidates, explicit = _discover(source, target)
    with ExitStack() as stack:
        result, notices = [], []
        for run in candidates:
            try:
                _path(run, manifest.MANIFEST)
                lock_path = run/'.run.lock'
                if lock_path.exists():
                    if lock_path.is_symlink():
                        raise ValueError('Source lock is a symlink')
                    lock = stack.enter_context(lock_path.open('rb'))
                    try:
                        fcntl.flock(lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
                    except BlockingIOError:
                        raise ValueError(f'Reuse source is active/locked: {run}') from None
                before = manifest.sha256(run/manifest.MANIFEST)
                value = manifest.load(run)
                manifest.validate_metadata(run, value)
                if manifest.sha256(run/manifest.MANIFEST) != before:
                    raise ValueError(f'Source manifest changed while reading: {run}')
                if not value.get('stage_records'):
                    notices.append(f'Source has no stage completion records; recompute or use its explicit legacy migration tool: {run}')
                result.append((run, value, before))
            except ERRORS as error:
                if explicit:
                    raise ValueError(f'Cannot read reuse source: {error}') from error
                notices.append(str(error))
        yield result, notices
        for run, value, checksum in result:
            if manifest.sha256(run/manifest.MANIFEST) != checksum:
                raise ValueError(f'Source manifest changed during reuse: {run}')


def _plan(value, sources, policy, target=None, notices=()):
    if policy not in POLICIES:
        raise ValueError(f'Unsupported reuse policy: {policy}')
    rows = []
    organisms = {o['slot']: o for o in value['organisms']}
    for kind, slots in (('gc', SLOTS), ('coding', SLOTS), ('train', TRAINS), ('pair', PAIRS)):
        for slot in slots:
            row = dict(kind=kind, slot=slot, action='recompute', reason='No compatible completed stage')
            if kind == 'coding' and organisms[slot]['gff_sha256'] is None:
                row.update(action='skip', reason='No GFF annotation')
                rows.append(row)
                continue
            if target is not None and _ready(target, value, kind, slot):
                row.update(action='skip', reason='Target stage already complete')
                rows.append(row)
                continue
            matches, rejected = [], []
            for run, old, checksum in sources:
                for key, record in old.get('stage_records', {}).items():
                    if record.get('kind') != kind:
                        continue
                    try:
                        source_slot = record['slot']
                        if key != f'{kind}:{source_slot}':
                            raise ValueError('Stage key disagrees with record slot')
                        # Validate against the original run before comparing logical target inputs.
                        allow_old = kind == 'pair' and policy == 'allow-existing-pairs'
                        _verified(run, old, kind, source_slot, allow_old)
                        if not _compatible(value, kind, slot, record, allow_old):
                            rejected.append('Stage inputs, options, implementation or LAST version differ')
                            continue
                        matches.append((run, record, checksum, source_slot))
                    except ERRORS as error:
                        rejected.append(str(error))
            # Equivalent logical slots within one run are a single source candidate.
            per_run = {}
            for match in sorted(matches, key=lambda m: (m[3] != slot, m[3])):
                per_run.setdefault(match[0], match)
            matches = list(per_run.values())
            if len(matches) > 1:
                raise ValueError(f'Ambiguous reusable stage {kind}:{slot}: '+', '.join(str(m[0]) for m in matches))
            if matches:
                run, record, checksum, source_slot = matches[0]
                row.update(action='reuse', reason='Verified matching completion record', source_run=str(run),
                           source_manifest_sha256=checksum, source_slot=source_slot, record=copy.deepcopy(record))
            elif rejected or notices:
                row['reason'] = '; '.join(dict.fromkeys(rejected or notices))
            rows.append(row)
    # A reused alignment carries the exact training profile that produced it.
    for row in rows:
        if row['kind'] != 'pair' or row['action'] != 'reuse':
            continue
        train = next(r for r in rows if r['kind'] == 'train' and r['slot'] == row['slot'])
        source = next(s for s in sources if str(s[0]) == row['source_run'])
        training = source[1]['stage_records'][f'train:{row["source_slot"]}']
        if train['action'] == 'skip' and target is not None:
            existing = value.get('stage_records', {}).get(f'train:{row["slot"]}')
            if existing and _file(existing, 'train')['sha256'] != _file(training, 'train')['sha256']:
                row.update(action='recompute', reason='Target training differs from pair training profile')
                continue
        train.update(action='reuse', reason='Exact training profile bundled with reused pair',
                     source_run=row['source_run'], source_manifest_sha256=row['source_manifest_sha256'],
                     source_slot=row['source_slot'], record=copy.deepcopy(training), bundle=True)
    return dict(policy=policy, stages=rows)


def plan(value, source=None, policy='same', target_run=None):
    if source is None:
        return _plan(value, [], policy, target_run)
    with _sources(source, target_run) as (sources, notices):
        return _plan(value, sources, policy, target_run, notices)


def _patterns():
    # These are strings in the repository's YAML, not shell commands.
    path = ROOT/'config/sbst_config.yaml'
    try:
        import yaml
        return yaml.safe_load(path.read_text())['patterns']
    except ImportError:
        result = {}
        for key in ('gc_content', 'coding_region_percentage', 'train', 'one2one'):
            result[key] = subprocess.check_output(['yq', 'eval', '-r', '.patterns.'+key, str(path)], text=True).strip()
        return result


def _destination(value, kind, slot, role, patterns):
    organisms = {o['slot']: o for o in value['organisms']}
    fields = {'date': value['date']}
    if kind in ('gc', 'coding'):
        fields['org_short'] = organisms[slot]['short_name']
        pattern = patterns['gc_content' if kind == 'gc' else 'coding_region_percentage']
        return 'statistics/misc/'+pattern.format(**fields)
    ref, query = slot.split('-')
    fields.update(org1_short=organisms[ref]['short_name'], org2_short=organisms[query]['short_name'])
    pattern = patterns['train' if role == 'train' else 'one2one']
    name = pattern.format(**fields)
    if role == 'maf_sorted':
        name = name[:-4]+'_sorted.maf' if name.endswith('.maf') else name+'_sorted.maf'
    return 'intermediateFiles/'+name


def apply(run, source, policy='same'):
    run = Path(run).resolve()
    _path(run, manifest.MANIFEST)
    value = manifest.load(run)
    patterns = _patterns()
    with _sources(source, run) as (sources, notices):
        result = _plan(value, sources, policy, run, notices)
        jobs, imported = {}, []
        for row in result['stages']:
            if row['action'] != 'reuse':
                continue
            record = copy.deepcopy(row['record'])
            record['slot'] = row['slot']
            record['provenance'] = dict(source_run=row['source_run'],
                source_manifest_sha256=row['source_manifest_sha256'], policy=policy,
                original_provenance=copy.deepcopy(record.get('provenance')))
            for output in record['outputs']:
                origin = _path(row['source_run'], output['path'])
                relative = _destination(value, row['kind'], row['slot'], output['role'], patterns)
                dest = _path(run, relative)
                expected = output['sha256']
                if dest.exists() and (not dest.is_file() or manifest.sha256(dest) != expected):
                    raise ValueError(f'Refusing to overwrite different target output: {dest}')
                if dest in jobs and jobs[dest][1] != expected:
                    raise ValueError(f'Conflicting imported target output: {dest}')
                jobs[dest] = (origin, expected)
                output['path'] = relative
            imported.append(record)
        # Stage every copy and verify every byte before publishing any output.
        staged = []
        try:
            for dest, (origin, expected) in jobs.items():
                if dest.exists():
                    continue
                dest.parent.mkdir(parents=True, exist_ok=True)
                fd, name = tempfile.mkstemp(prefix=dest.name+'.reuse.', dir=dest.parent)
                os.close(fd)
                temporary = Path(name)
                staged.append((temporary, dest))
                shutil.copyfile(origin, temporary)
                if manifest.sha256(temporary) != expected or manifest.sha256(origin) != expected:
                    raise ValueError(f'Source changed during copy: {origin}')
            for original, old, checksum in sources:
                if manifest.sha256(original/manifest.MANIFEST) != checksum:
                    raise ValueError(f'Source manifest changed during reuse: {original}')
            for temporary, dest in staged:
                # Link a newly copied temporary inode: atomic and cannot overwrite a racer.
                try:
                    os.link(temporary, dest)
                except FileExistsError:
                    if manifest.sha256(dest) != jobs[dest][1]:
                        raise ValueError(f'Refusing to overwrite changed target: {dest}') from None
            for record in imported:
                value.setdefault('stage_records', {})[f'{record["kind"]}:{record["slot"]}'] = record
                for output in record['outputs']:
                    if output['role'] == 'train':
                        manifest.seal_file(_path(run, output['path']))
                _record_artifacts(value, run, record)
            if imported:
                manifest.atomic_json(run/manifest.MANIFEST, value)
        finally:
            for temporary, dest in staged:
                temporary.unlink(missing_ok=True)
    return result


def relocate(run, old, new):
    """Record a successful maf-sort relocation, retaining the producing provenance."""
    run = Path(run).resolve()
    _path(run, manifest.MANIFEST)
    value = manifest.load(run)
    destination = _path(run, new)
    validate_maf(destination)
    changed = False
    for record in value.get('stage_records', {}).values():
        if record.get('kind') != 'pair':
            continue
        for output in record['outputs']:
            if output['path'] == old and output['role'].startswith('maf'):
                output.update(role='maf_sorted', path=new, sha256=manifest.sha256(destination))
                _record_artifacts(value, run, record)
                changed = True
    if changed:
        value['artifacts'] = [a for a in value['artifacts'] if a['path'] != old]
        manifest.atomic_json(run/manifest.MANIFEST, value)
    return changed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('register'); p.add_argument('run'); p.add_argument('kind'); p.add_argument('slot'); p.add_argument('outputs', nargs='+')
    p = sub.add_parser('ready'); p.add_argument('run'); p.add_argument('kind'); p.add_argument('slot')
    p = sub.add_parser('import'); p.add_argument('run'); p.add_argument('source'); p.add_argument('--policy', choices=POLICIES, default='same')
    p = sub.add_parser('plan'); p.add_argument('--descriptor', required=True); p.add_argument('source'); p.add_argument('--policy', choices=POLICIES, default='same')
    p = sub.add_parser('relocate'); p.add_argument('run'); p.add_argument('old'); p.add_argument('new')
    args = parser.parse_args()
    try:
        if args.command == 'ready':
            return 0 if ready(args.run, args.kind, args.slot) else 1
        if args.command == 'register':
            pairs = [o.split('=', 1) for o in args.outputs]
            if any(len(o) != 2 for o in pairs) or len(dict(pairs)) != len(pairs):
                raise ValueError('Outputs must be unique ROLE=PATH values')
            result = register(args.run, args.kind, args.slot, dict(pairs))
        elif args.command == 'import':
            result = apply(args.run, args.source, args.policy)
        elif args.command == 'plan':
            result = plan(json.loads(args.descriptor), args.source, args.policy)
        else:
            result = relocate(args.run, args.old, args.new)
        if args.command=='register':
            print('Completed '+args.kind+':'+args.slot)
        elif args.command=='relocate':
            if result:print('Recorded sorted alignment: '+args.new)
        else:
            shown=copy.deepcopy(result)
            for row in shown['stages']:
                record=row.pop('record',None)
                if record:row['producer_last_version']=record['producer'].get('last_version')
            print(json.dumps(shown,sort_keys=True))
        return 0
    except ERRORS as error:
        print(f'Stage reuse error: {error}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
