#!/usr/bin/env python3
"""Versioned, content-addressed run provenance; filenames are never identities."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import platform
import re
from datetime import datetime, timezone

SCHEMA = 2
ALGORITHM = 'pipeline-v2-canonical-subratio'
MANIFEST = 'metadata/metadata_manifest.json'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f'.tmp-{os.getpid()}')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')
    temporary.replace(path)


def load(run):
    path = Path(run) / MANIFEST
    if not path.is_file():
        raise ValueError(f'New-format manifest required; rerun the analysis: {run}')
    value = json.loads(path.read_text())
    if value.get('schema_version') != SCHEMA:
        raise ValueError(f'Unsupported legacy manifest; rerun the analysis: {path}')
    required = ('run_key', 'trio_key', 'organisms', 'artifacts', 'status', 'created_at', 'options')
    if not all(key in value for key in required):
        raise ValueError(f'Incomplete v2 manifest: {path}')
    slots = [entry.get('slot') for entry in value['organisms']]
    if sorted(slots) != ['org1', 'org2', 'org3']:
        raise ValueError(f'Invalid organism slots: {path}')
    if value['status'] not in ('pending', 'identity_complete', 'counts_complete', 'complete'):
        raise ValueError(f'Invalid run status: {path}')
    paths = [a['path'] for a in value['artifacts']]
    if len(paths) != len(set(paths)):
        raise ValueError(f'Duplicate artifact entries: {path}')
    logical = [(a['kind'],a['slot']) for a in value['artifacts'] if a['kind'] != 'auxiliary']
    if len(logical) != len(set(logical)):
        raise ValueError(f'Duplicate logical artifact: {path}')
    for entry in value['organisms']:
        if not entry.get('fasta_sha256') or 'gff_sha256' not in entry:
            raise ValueError(f'Incomplete input provenance: {path}')
    return value


def validated(run, complete=False):
    value = load(run)
    validate_contents(run, value, complete)
    return value


def validate_contents(run, value, complete=False):
    if complete and value['status'] != 'complete':
        raise ValueError(f'Completed analysis required: {run}')
    if complete:
        required = {('train', pair) for pair in ('org1-org2','org1-org3','org2-org3')}
        required |= {('gc', slot) for slot in ('org1','org2','org3')}
        required.add(('sbst_ratio','all'))
        outgroup = next(e for e in value['organisms'] if e['slot']=='org1')
        kinds = ['single','dinuc']
        if outgroup['gff_sha256'] is not None:
            kinds += ['single_ncds','dinuc_ncds']
        required |= {(kind,slot) for kind in kinds for slot in ('org2','org3')}
        for entry in value['organisms']:
            if entry['gff_sha256'] is not None:
                required.add(('coding',entry['slot']))
        present = {(a['kind'],a['slot']) for a in value['artifacts']}
        if not required.issubset(present):
            raise ValueError(f'Missing completed artifacts: {sorted(required-present)}')
        validate_plots(run, value)
    validate_metadata(run,value)
    for entry in value['artifacts']:
        path = artifact_path(run, entry['path'])
        if not path.is_file() or sha256(path) != entry['sha256']:
            raise ValueError(f'Artifact checksum mismatch: {path}')
        if entry['kind'] == 'train':
            from train_identity import parse_identity
            measured = parse_identity(path)
            if measured != entry.get('identity_percent'):
                raise ValueError(f'Training identity disagrees with provenance: {path}')


def validate_metadata(run,value):
    selection=value.get('selection_provenance')
    if selection and sha256(artifact_path(run,selection['path'])) != selection['sha256']:
        raise ValueError('Selection provenance checksum mismatch')
    for entry in value['organisms']:
        for field in ('metadata_json','taxonomy_json'):
            if entry.get(field):
                path=artifact_path(run,entry[field])
                if sha256(path) != entry.get(field+'_sha256'):
                    raise ValueError(f'Metadata checksum mismatch: {path}')


def validate_plots(run, value):
    records = {(a['kind'], a['slot']): a for a in value['artifacts']}
    for slot in ('org2', 'org3'):
        for noncoding in (False, True):
            channel = 'single_ncds' if noncoding else 'single'
            if (channel, slot) not in records:
                if noncoding:
                    continue
                raise ValueError(f'Missing single-base count artifact: {slot}')
            prefix = 'pdf_ncds_' if noncoding else 'pdf_'
            required = [prefix + kind for kind in ('norm', 'sbst', 'ori', 'logratio')]
            required.append('pdf_dinuc_ncds_tsv' if noncoding else 'pdf_dinuc_tsv')
            for kind in required:
                if (kind, slot) not in records:
                    raise ValueError(f'Missing completed plot: {slot} {kind}')
                path = artifact_path(run, records[(kind, slot)]['path'])
                data = path.read_bytes()
                if not data.startswith(b'%PDF-') or b'%%EOF' not in data[-1024:]:
                    raise ValueError(f'Incomplete PDF: {path}')


def discover(root):
    root = Path(root).resolve()
    candidates = []
    visited = set()
    for current, directories, files in os.walk(root, followlinks=True):
        actual = Path(current).resolve()
        if actual in visited:
            directories[:] = []
            continue
        visited.add(actual)
        if 'metadata_manifest.json' in files and Path(current).name == 'metadata':
            candidates.append(Path(current).parent.resolve())
    candidates = sorted(set(candidates))
    if (root / MANIFEST).is_file():
        candidates = [root]
    if not candidates:
        raise ValueError(f'No v2 runs found; rerun legacy analyses: {root}')
    latest = {}
    for run in candidates:
        value = load(run)
        if value['status'] != 'complete':
            continue
        key = value['trio_key']
        timestamp = value['created_at']
        if key not in latest or timestamp > latest[key][0]:
            latest[key] = (timestamp, run)
    return [item[1] for _, item in sorted(latest.items())]


def artifacts(run, kind, slot=None, validated_manifest=None):
    value = validated_manifest if validated_manifest is not None else validated(run)
    return [artifact_path(run, a['path']) for a in value['artifacts']
            if a['kind'] == kind and (slot is None or a['slot'] == slot)]


def artifact_path(run, relative):
    path = (Path(run) / relative).resolve()
    if not path.is_relative_to(Path(run).resolve()):
        raise ValueError(f'Artifact path escapes run: {relative}')
    return path


def effective_last_options():
    from last_options import config
    return config()


def descriptor(args):
    if not re.fullmatch(r'[A-Za-z0-9_.-]+',args.label) or args.label in ('.','..'):
        raise ValueError('Run label must be a single filename component')
    organisms = []
    extra_path = os.environ.get('EVOSUBSTER_METADATA_INPUT')
    extra = json.loads(Path(extra_path).read_text()) if extra_path else [{}] * 3
    for i in range(3):
        fasta, gff, accession, name = args.inputs[i * 4:(i + 1) * 4]
        entry = dict(extra[i])
        entry.update(slot=f'org{i+1}', role='outgroup' if i == 0 else 'ingroup',
                     accession=accession or None, short_name=name,
                     fasta_path=str(Path(fasta).resolve()), fasta_sha256=sha256(fasta),
                     gff_path=None if gff == 'NO_GFF_FILE' else str(Path(gff).resolve()),
                     gff_sha256=None if gff == 'NO_GFF_FILE' else sha256(gff))
        organisms.append(entry)
    try:
        version = subprocess.check_output(['lastal', '--version'], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        raise ValueError('Cannot identify LAST version')
    identity = [{k: entry[k] for k in ('slot', 'accession', 'fasta_sha256', 'gff_sha256')}
                for entry in organisms]
    config_root = Path(__file__).resolve().parents[2] / 'config'
    options = dict(last=effective_last_options(),
                   python_version=platform.python_version(),
                   r_version=subprocess.check_output(['Rscript','--version'],text=True,stderr=subprocess.STDOUT).strip(),
                   threads=args.threads, last_version=version, algorithm=ALGORITHM,
                   config={p.name: sha256(p) for p in sorted(config_root.glob('*.yaml'))})
    trio_key = digest(identity)
    run_key = digest(dict(inputs=identity, options=options))
    selection = os.environ.get('EVOSUBSTER_SELECTION_SOURCE')
    provenance = None
    if selection and Path(selection).is_file():
        provenance = dict(source=str(Path(selection).resolve()),sha256=sha256(selection),path='metadata/selected_trios.tsv')
    packages = subprocess.check_output(['Rscript','-e',
        'p <- c("sysfonts","showtext","ggplot2"); cat(paste(p,vapply(p,function(x) if(requireNamespace(x,quietly=TRUE)) as.character(packageVersion(x)) else "missing",character(1)),sep="=",collapse=";"))'],text=True).strip()
    return dict(plot_package_versions=packages,selection_provenance=provenance,schema_version=SCHEMA, algorithm=ALGORITHM, trio_key=trio_key,
                run_key=run_key, date=args.label, options=options, organisms=organisms,
                artifacts=[], status='pending', created_at=datetime.now(timezone.utc).isoformat())


def recover_staged_outputs(run, previous, current):
    """Recover sealed stages written after the last successful pair checkpoint.

    A mere filename or MAF header never qualifies. The stage seal must match the
    current LAST executable, full settings, exact genomes and checkpointed train.
    Failed-stage stderr/timing survive separately from reusable analysis outputs.
    """
    if os.environ.get('LASTAL_PIPELINE_MODE') != 'staged':
        return
    from last_options import config, output_seal
    if previous['options'].get('last') != config():
        return
    organisms = {entry['slot']:entry for entry in current['organisms']}
    names = {entry['slot']:entry.get('short_name') for entry in previous['organisms']}
    if not all(names.values()) or not previous.get('date'):
        return
    existing = {a['path']:a for a in previous['artifacts']}
    for query in ('org2','org3'):
        prefix = f"intermediateFiles/{names['org1']}2{names[query]}"
        train = f"{prefix}_{previous['date']}.train"
        maf = f"{prefix}_one2one_{previous['date']}.maf"
        for stage in ('lastal','last-split','maf-linked'):
            for kind in ('stderr','time'):
                pattern = Path(maf).name + f'.{stage}.{kind}.*'
                for diagnostic in (run/'intermediateFiles').glob(pattern):
                    if not diagnostic.is_file() or diagnostic.is_symlink() or not diagnostic.name.rsplit('.',1)[-1].isdigit():
                        continue
                    checksum = sha256(diagnostic)
                    archive = run/'metadata/diagnostics'/diagnostic.name
                    if archive.exists() and sha256(archive) != checksum:
                        archive = archive.with_name(archive.name+'.'+checksum[:16])
                    archive.parent.mkdir(parents=True, exist_ok=True)
                    old_name = diagnostic.relative_to(run).as_posix()
                    diagnostic.replace(archive)
                    existing.pop(old_name,None)
                    relative = archive.relative_to(run).as_posix()
                    existing[relative] = dict(kind='auxiliary',slot=None,path=relative,sha256=checksum)
        training = existing.get(train)
        if not training or not (run/train).is_file() or sha256(run/train) != training['sha256']:
            continue
        inputs = [organisms['org1'].get('fasta_path'),organisms[query].get('fasta_path'),str(run/train)]
        if not all(inputs):
            continue
        for relative in (maf+'.lastal.maf',maf+'.raw',maf):
            seal = relative+'.complete.json'
            if not (run/seal).is_file():
                continue
            try:
                output_seal(run/relative,inputs,verify=True)
            except (ValueError,OSError,KeyError,TypeError):
                existing.pop(relative,None)
                existing.pop(seal,None)
                continue
            for filename in (relative,seal):
                existing[filename] = dict(kind='auxiliary',slot=None,path=filename,sha256=sha256(run/filename))
            guard = maf+'.last_options.json'
            if (run/guard).is_file() and json.loads((run/guard).read_text()) == previous['options']['last']:
                existing[guard] = dict(kind='auxiliary',slot=None,path=guard,sha256=sha256(run/guard))
    previous['artifacts'] = list(existing.values())


def prepare(run, value):
    run = Path(run)
    path = run / MANIFEST
    if path.exists():
        previous = load(run)
        if previous['run_key'] != value['run_key']:
            raise ValueError('Existing run inputs/options disagree with manifest')
        validate_metadata(run, previous)
        recover_staged_outputs(run, previous, value)
        value = previous
        previous_count = len(previous['artifacts'])
        verified = {a['path'] for a in value['artifacts']
                    if artifact_path(run, a['path']).is_file()
                    and sha256(artifact_path(run, a['path'])) == a['sha256']}
        # Only pipeline-owned generated trees are considered for recovery.
        for directory in ('intermediateFiles', 'statistics', 'figs'):
            for item in (run / directory).rglob('*'):
                if item.is_file() or item.is_symlink():
                    if item.relative_to(run).as_posix() not in verified:
                        item.unlink()
            for item in sorted((run / directory).rglob('*'), key=lambda p: len(p.parts), reverse=True):
                if item.is_dir() and not any(item.iterdir()):
                    item.rmdir()
        value['artifacts'] = [a for a in value['artifacts'] if a['path'] in verified]
        if len(verified) != previous_count:
            value['status'] = 'pending'
    elif run.exists() and any(p.name != '.run.lock' for p in run.iterdir()):
        raise ValueError('Existing output has no v2 manifest; use a new output directory')
    for entry in value['organisms']:
        if not entry.get('metadata_json'):
            target = run / 'metadata' / f"{entry['slot']}_assembly.json"
            atomic_json(target, {'local_input':True,'display_name':entry.get('short_name',entry['slot'])})
            entry['metadata_json'] = str(target.resolve())
        for field in ('metadata_json', 'taxonomy_json'):
            source = entry.get(field)
            if source and not Path(source).is_absolute() and path.exists():
                source = str(run / source)
            if source and Path(source).is_file():
                suffix = 'taxonomy' if field == 'taxonomy_json' else 'assembly'
                target = run / 'metadata' / f"{entry['slot']}_{suffix}.json"
                target.parent.mkdir(parents=True, exist_ok=True)
                if Path(source).resolve() != target.resolve():
                    shutil.copyfile(source, target)
                entry[field] = target.relative_to(run).as_posix()
                entry[field + '_sha256'] = sha256(target)
    selection = value.get('selection_provenance')
    if selection and not (run / selection['path']).exists():
        source=Path(selection['source'])
        if sha256(source) != selection['sha256']:
            raise ValueError('Selection provenance changed before publication')
        shutil.copyfile(source,run / selection['path'])
    atomic_json(path, value)


def checkpoint(run, status, specs):
    run = Path(run)
    value = load(run)
    # Successful later stages may consume their sealed intermediates.
    existing = {a['path']: a for a in value['artifacts']
                if artifact_path(run, a['path']).is_file()}
    for spec in specs:
        kind, slot, filename = spec.split(':', 2)
        path = artifact_path(run, filename)
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f'Required output missing/empty: {path}')
        artifact = dict(kind=kind, slot=slot, path=filename, sha256=sha256(path))
        if kind == 'train':
            from train_identity import parse_identity
            artifact['identity_percent'] = parse_identity(path)
            artifact['reference_slot'], artifact['query_slot'] = slot.split('-')
        existing[filename] = artifact
        if kind == 'sbst_ratio' or (kind in ('single','single_ncds','dinuc','dinuc_ncds') and slot == 'org2'):
            marker = Path(str(path) + '.complete')
            if marker.is_file():
                relative = marker.relative_to(run).as_posix()
                existing[relative] = dict(kind='auxiliary',slot=None,path=relative,sha256=sha256(marker))
    # Snapshot all intermediate files only after the stage succeeds.
    for directory in ('intermediateFiles',):
        for path in (run / directory).rglob('*'):
            if path.is_file() and not path.is_symlink() and not any(
                token in path.name for token in ('.tmp', '.lock', '.cache-copy')
            ):
                relative = path.relative_to(run).as_posix()
                if relative not in existing:
                    existing[relative] = dict(kind='auxiliary', slot=None, path=relative, sha256=sha256(path))
    if status == 'complete':
        plots = json.loads((run / 'metadata/plot_completion.json').read_text())
        if plots.get('status') != 'complete' or not plots.get('artifacts'):
            raise ValueError('Complete plot inventory required')
        for entry in plots['artifacts']:
            path = artifact_path(run, entry['path'])
            if sha256(path) != entry['sha256']:
                raise ValueError(f'Plot checksum mismatch: {path}')
            existing[entry['path']] = entry
    value['artifacts'] = list(existing.values())
    order = ['pending', 'identity_complete', 'counts_complete', 'complete']
    value['status'] = max((value['status'], status), key=order.index)
    value['updated_at'] = datetime.now(timezone.utc).isoformat()
    if status == 'complete':
        validate_contents(run, value, complete=True)
    atomic_json(run / MANIFEST, value)


def training_directory(base, ref, query, ref_acc, query_acc):
    value = dict(algorithm='training-v2', reference=sha256(ref), query=sha256(query),
                 reference_accession=ref_acc, query_accession=query_acc,
                 last=effective_last_options(),
                 threads=os.environ.get('THREAD_NUM', '8'),
                 last_version=subprocess.check_output(['lastal', '--version'], text=True).strip())
    return Path(base).resolve() / 'v2' / digest(value)


def bundle_save(directory, accession, fasta):
    directory = Path(directory).resolve()
    fasta = Path(fasta).resolve()
    gff = directory / 'genomic.gff'
    if gff.is_file():
        with fasta.open() as handle:
            sequence_ids = {line[1:].split()[0] for line in handle if line.startswith('>')}
        with gff.open() as handle:
            annotated = {line.split('\t', 1)[0] for line in handle
                         if line.strip() and not line.startswith('#')}
        if not annotated.issubset(sequence_ids):
            raise ValueError('GFF sequence IDs disagree with downloaded FASTA')
    atomic_json(directory / 'assembly_bundle.json', dict(schema_version=2, accession=accession,
                fasta=fasta.name, fasta_sha256=sha256(fasta),
                gff=gff.name if gff.is_file() else None,
                gff_sha256=sha256(gff) if gff.is_file() else None))


def bundle_gff(fasta, expected_accession=None):
    fasta = Path(fasta).resolve()
    path = fasta.parent / 'assembly_bundle.json'
    if not path.is_file():
        raise ValueError(f'Verified assembly bundle required: {path}')
    value = json.loads(path.read_text())
    if expected_accession and value.get('accession') != expected_accession:
        raise ValueError('Assembly accession disagrees with bundle')
    if value.get('schema_version') != 2 or value['fasta'] != fasta.name or value['fasta_sha256'] != sha256(fasta):
        raise ValueError('FASTA disagrees with assembly bundle')
    if value['gff'] is None:
        return 'NO_GFF_FILE'
    gff = fasta.parent / value['gff']
    if sha256(gff) != value['gff_sha256']:
        raise ValueError('GFF disagrees with assembly bundle')
    return str(gff)


def seal_directory(path, verify=False):
    path = Path(path)
    seal = path / '.database-complete.json'
    files = {p.relative_to(path).as_posix(): sha256(p) for p in path.rglob('*')
             if p.is_file() and p != seal and not p.name.endswith('.lock')}
    if not files:
        raise ValueError(f'Empty database: {path}')
    if verify:
        if json.loads(seal.read_text()) != files:
            raise ValueError(f'Incomplete or modified database: {path}')
    else:
        atomic_json(seal, files)


def seal_file(path, verify=False):
    path = Path(path)
    from train_identity import parse_identity
    parse_identity(path)
    value = dict(sha256=sha256(path))
    seal = Path(str(path) + '.complete.json')
    if verify:
        if json.loads(seal.read_text()) != value:
            raise ValueError(f'Modified training cache: {path}')
    else:
        atomic_json(seal, value)


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('path', 'prepare'):
        p = sub.add_parser(name)
        p.add_argument('base'); p.add_argument('label'); p.add_argument('threads', type=int)
        p.add_argument('inputs', nargs=12)
    p = sub.add_parser('bundle-save'); p.add_argument('values', nargs=3)
    p = sub.add_parser('bundle-gff'); p.add_argument('fasta'); p.add_argument('accession', nargs='?')
    p = sub.add_parser('discover'); p.add_argument('root')
    p = sub.add_parser('read'); p.add_argument('run')
    for name in ('seal-directory', 'verify-directory', 'seal-file', 'verify-file'):
        p = sub.add_parser(name); p.add_argument('path')
    p = sub.add_parser('training-directory'); p.add_argument('values', nargs=5)
    p = sub.add_parser('checkpoint'); p.add_argument('run'); p.add_argument('status')
    p.add_argument('specs', nargs='*')
    args = parser.parse_args()
    if args.command in ('seal-directory', 'verify-directory'):
        seal_directory(args.path, args.command.startswith('verify'))
    elif args.command in ('seal-file', 'verify-file'):
        seal_file(args.path, args.command.startswith('verify'))
    elif args.command == 'bundle-save':
        bundle_save(*args.values)
    elif args.command == 'bundle-gff':
        print(bundle_gff(args.fasta, args.accession))
    elif args.command == 'read':
        print(json.dumps(validated(args.run, complete=True)))
    elif args.command == 'discover':
        print(json.dumps([str(p) for p in discover(args.root)]))
    elif args.command == 'training-directory':
        print(training_directory(*args.values))
    elif args.command == 'checkpoint':
        checkpoint(args.run, args.status, args.specs)
    else:
        value = descriptor(args)
        run = Path(args.base).resolve() / value['trio_key'] / f"{args.label}-{value['run_key']}"
        if args.command == 'prepare':
            prepare(run, value)
        print(run)


if __name__ == '__main__':
    main()
