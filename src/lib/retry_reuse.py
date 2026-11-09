#!/usr/bin/env python3
"""Explicit, copy-only import of certified LAST 1648 outputs into a 1654 retry.

Legacy alignments were not checkpointed until counting finished. Their certification
therefore also requires the exact successful commands in the original run's log.
This is deliberately separate from the normal settings and completion guards.
"""
import argparse
import copy
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import sys

import run_manifest as manifest
from train_identity import parse_identity

PAIRS = ('org1-org2', 'org1-org3')
ALIGNMENT_TRIOS = {
    ('GCF_040869185.1', 'GCA_026956025.1', 'GCF_049350105.2'),
    ('GCF_964188355.1', 'GCF_900496995.4', 'GCF_026419915.1'),
    ('GCA_910950805.1', 'GCF_001522545.3', 'GCA_038429875.2'),
}
POLICY = 'legacy-1648-to-1654-certified-pairs-v1'


def _without_threads(words):
    result = []
    index = 0
    while index < len(words):
        word = words[index]
        if word == '-P':
            if index + 1 >= len(words):
                raise ValueError('Missing LAST thread value')
            index += 2
        elif re.fullmatch(r'-P[0-9]+', word):
            index += 1
        else:
            result.append(word)
            index += 1
    return result


def compatible(target, source):
    fields = ('slot', 'accession', 'fasta_sha256', 'gff_sha256', 'short_name')
    identity = lambda value: [{k: o.get(k) for k in fields}
                              for o in sorted(value['organisms'], key=lambda o:o['slot'])]
    if identity(target) != identity(source) or target.get('date') != source.get('date'):
        raise ValueError('Reuse inputs, accessions, display names or date disagree')
    old, new = source['options'], target['options']
    if old.get('last_version') != 'lastal 1648' or new.get('last_version') != 'lastal 1654':
        raise ValueError('Reuse policy only permits LAST 1648 to 1654')
    if any(old.get(k) != new.get(k) for k in ('algorithm', 'config')):
        raise ValueError('Reuse algorithm/configuration disagree')
    oldlast, newlast = old.get('last', {}), new.get('last', {})
    keys = {'schema', 'lastdb_train', 'lastdb_align', 'last_train', 'lastal', 'last_split'}
    if set(oldlast) != keys or set(newlast) != keys or oldlast['schema'] != newlast['schema']:
        raise ValueError('Unrecognized LAST configuration schema')
    if any(_without_threads(oldlast[k]) != _without_threads(newlast[k]) for k in keys - {'schema'}):
        raise ValueError('Reuse LAST semantic parameters disagree')


def validate_maf(path):
    """Stream a nonempty pairwise MAF; reject truncated rows and invalid coordinates."""
    blocks, rows, active = 0, [], False

    def finish():
        if not active:
            return
        if len(rows) != 2 or len(rows[0]) != len(rows[1]):
            raise ValueError(f'Incomplete pairwise MAF block: {path}')

    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split()
            if fields[0] == 'a':
                finish()
                blocks += 1
                active, rows = True, []
            elif fields[0] == 's':
                if not active or len(fields) != 7:
                    raise ValueError(f'Invalid MAF sequence row {path}:{number}')
                _, name, start, size, strand, source_size, sequence = fields
                try:
                    start, size, source_size = int(start), int(size), int(source_size)
                except ValueError:
                    raise ValueError(f'Invalid MAF coordinate {path}:{number}') from None
                if (start < 0 or size <= 0 or source_size < 1 or start + size > source_size
                        or strand not in ('+', '-') or len(sequence.replace('-', '')) != size
                        or not re.fullmatch('[A-Za-z*-]+', sequence)):
                    raise ValueError(f'Invalid MAF coordinate/length {path}:{number}')
                rows.append(sequence)
                if len(rows) > 2:
                    raise ValueError(f'MAF is not pairwise: {path}')
            elif fields[0] not in ('p', 'q', 'i') or not active:
                raise ValueError(f'Unsupported MAF row {path}:{number}')
    finish()
    if not blocks:
        raise ValueError(f'Empty MAF: {path}')
    return blocks


def _tokens(command):
    lex = shlex.shlex(command, posix=True, punctuation_chars='|>')
    lex.whitespace_split = True
    words = list(lex)
    normalized = []
    for word in words:
        if re.fullmatch(r'-P[0-9]+', word):
            normalized.extend(('-P', word[2:]))
        elif word.startswith('--') and '=' in word:
            normalized.extend(word.split('=', 1))
        else:
            normalized.append(word)
    return normalized


def _completed_commands(log, source):
    """Only commands in this run's pwd section with all three timing rows qualify."""
    commands, active = [], False
    lines = Path(log).read_text().splitlines()
    for index, line in enumerate(lines):
        if line.startswith('pwd: '):
            active = Path(line[5:]).resolve() == Path(source).resolve()
        if not active or not line.startswith('time '):
            continue
        trailing = []
        for next_line in lines[index+1:]:
            if next_line.startswith(('time ', 'bash ', '---', 'pwd: ', '===', 'cleanup:')):
                break
            trailing.append(next_line)
        evidence = '\n'.join(trailing)
        if re.search(r'killed|terminated|segmentation|skipping sequence|error:|failed|exit=[1-9]', evidence, re.I):
            continue
        if all(re.search(r'^'+name+r'\s+[0-9]+m[0-9]+(?:\.[0-9]+)?s\s*$', evidence, re.M)
               for name in ('real', 'user', 'sys')):
            commands.append(dict(line=index+1, command=line, argv=_tokens(line), timing=evidence))
    return commands


def _pair_paths(value, pair):
    organisms = {o['slot']:o for o in value['organisms']}
    ref, query = (organisms[s] for s in pair.split('-'))
    prefix = f"intermediateFiles/{ref['short_name']}2{query['short_name']}"
    return f"{prefix}_{value['date']}.train", f"{prefix}_one2one_{value['date']}.maf"


def _pair_evidence(value, pair, commands, sorted_maf=False):
    train, maf = _pair_paths(value, pair)
    organisms = {o['slot']:o for o in value['organisms']}
    ref, query = (organisms[s] for s in pair.split('-'))
    settings = value['options']['last']
    db_args = settings['lastdb_align']
    seed = db_args[db_args.index('-u')+1]
    if not re.fullmatch('[A-Za-z0-9_.-]+', seed):
        raise ValueError('Legacy custom seed files cannot be imported')
    db = f"{ref['short_name']}db_{value['date']}.seed-{seed}"
    prefix = ['time','lastal','-p',train] + settings['lastal'] + [f'intermediateFiles/{db}/{db}']
    suffix = ['|','last-split'] + settings['last_split'] + ['>',maf+'.raw']
    align = [c for c in commands if c['argv'][:len(prefix)] == prefix
             and c['argv'][len(prefix)+1:] == suffix
             and Path(c['argv'][len(prefix)]).name == Path(query['fasta_path']).name
             and re.match(re.escape(query['accession'])+r'(?:_|\.)', Path(c['argv'][len(prefix)]).name)]
    linked = [c for c in commands if c['argv'] == ['time','maf-linked',maf+'.raw','>',maf]]
    if len(align) != 1 or len(linked) != 1 or linked[0]['line'] <= align[0]['line']:
        raise ValueError(f'Missing unambiguous completed alignment commands for {pair}')
    result = [align[0], linked[0]]
    if sorted_maf:
        sorted_path = maf[:-4]+'_sorted.maf'
        sort = [c for c in commands if c['argv'] == ['time','maf-sort',maf,'>',sorted_path]]
        if len(sort) != 1 or sort[0]['line'] <= linked[0]['line']:
            raise ValueError(f'Missing completed maf-sort for {pair}')
        result.append(sort[0])
    return [{k:v for k,v in c.items() if k != 'argv'} for c in result]


def _verified_entry(run, entry):
    path = manifest.artifact_path(run, entry['path'])
    if not path.is_file() or (Path(run)/entry['path']).is_symlink() or manifest.sha256(path) != entry['sha256']:
        raise ValueError(f'Source checksum mismatch: {entry["path"]}')
    return path


def _copy(target, source, relative, expected):
    original = manifest.artifact_path(source, relative)
    dest = manifest.artifact_path(target, relative)
    if dest.exists():
        if manifest.sha256(dest) != expected:
            raise ValueError(f'Refusing to overwrite existing retry output: {dest}')
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    temporary = dest.with_name(dest.name+f'.reuse.tmp.{os.getpid()}')
    try:
        shutil.copyfile(original, temporary)
        if manifest.sha256(temporary) != expected:
            raise ValueError(f'Source changed during copy: {original}')
        temporary.replace(dest)
    finally:
        if temporary.exists(): temporary.unlink()


def pair_ready(target, pair):
    """Fail closed; no imported profile/alignment bypass without intact provenance."""
    try:
        value = manifest.load(target)
        reuse = value.get('legacy_reuse', {})
        if reuse.get('policy') != POLICY or reuse.get('target_run_key') != value['run_key']:
            return False
        record = reuse.get('pairs', {}).get(pair)
        if not record:
            return False
        current = {a['path']:a for a in value['artifacts']}
        for entry in record['files']:
            if current.get(entry['path'], {}).get('sha256') != entry['sha256']:
                return False
            _verified_entry(target, entry)
        manifest.seal_file(manifest.artifact_path(target, record['train']), verify=True)
        return True
    except (ValueError, OSError, KeyError, TypeError):
        return False


def current_pair_ready(target, pair):
    """Accept a native sorted pair only after its successful manifest checkpoint."""
    try:
        value = manifest.load(target)
        if pair in value.get('legacy_reuse', {}).get('pairs', {}):
            return False
        train, maf = _pair_paths(value, pair)
        current = {a['path']:a for a in value['artifacts']}
        names = (train, train+'.complete.json', maf[:-4]+'_sorted.maf', maf+'.last_options.json')
        for name in names:
            _verified_entry(target, current[name])
        if json.loads(manifest.artifact_path(target, names[-1]).read_text()) != value['options']['last']:
            return False
        manifest.seal_file(manifest.artifact_path(target, train), verify=True)
        return True
    except (ValueError, OSError, KeyError, TypeError):
        return False


def import_run(target, source, log):
    target, source, log = Path(target).resolve(), Path(source).resolve(), Path(log).resolve()
    if target == source or target.is_relative_to(source) or source.is_relative_to(target):
        raise ValueError('Reuse source and target must be separate run directories')
    new, old = manifest.load(target), manifest.load(source)
    compatible(new, old)
    manifest.validate_metadata(source, old)
    record = dict(policy=POLICY, source_run=str(source), source_run_key=old['run_key'],
                  source_manifest_sha256=manifest.sha256(source/manifest.MANIFEST),
                  source_log=str(log), source_log_sha256=manifest.sha256(log),
                  source_options=copy.deepcopy(old['options']), target_run_key=new['run_key'],
                  pairs={}, skipped=[])
    previous = new.get('legacy_reuse')
    if previous:
        if any(previous.get(k) != record[k] for k in ('policy','source_run','source_manifest_sha256',
                                                      'source_log','source_log_sha256','target_run_key')):
            raise ValueError('Retry source provenance changed')
        for pair in previous['pairs']:
            if not pair_ready(target,pair):
                raise ValueError(f'Imported pair checksum mismatch: {pair}')
        return previous
    existing = {a['path']:a for a in new['artifacts']}
    for entry in old['artifacts']:
        if entry['kind'] not in ('gc','coding'):
            continue
        try:
            _verified_entry(source,entry)
        except (OSError,ValueError) as error:
            record['skipped'].append(str(error)); continue
        if entry['path'] in existing:
            _verified_entry(target,existing[entry['path']]); continue
        _copy(target,source,entry['path'],entry['sha256'])
        existing[entry['path']] = copy.deepcopy(entry)
    commands = _completed_commands(log,source)
    accessions = tuple(o['accession'] for o in sorted(old['organisms'], key=lambda o:o['slot']))
    for pair in PAIRS if accessions in ALIGNMENT_TRIOS else ():
        train, maf = _pair_paths(old,pair)
        sorted_maf = maf[:-4]+'_sorted.maf'
        chosen = sorted_maf if (source/sorted_maf).is_file() else maf
        if not (source/chosen).is_file():
            continue
        # A retry's own completed pair must never be replaced by legacy output.
        if any(name in existing for name in (maf,sorted_maf)):
            continue
        try:
            training = next(a for a in old['artifacts'] if a['kind']=='train' and a['slot']==pair and a['path']==train)
            _verified_entry(source,training)
            if parse_identity(source/train) != training.get('identity_percent'):
                raise ValueError(f'Source identity disagrees: {pair}')
            manifest.seal_file(source/train,verify=True)
            evidence = _pair_evidence(old,pair,commands,chosen==sorted_maf)
            maf_path = manifest.artifact_path(source,chosen)
            if (source/chosen).is_symlink(): raise ValueError('Symlink MAF is not eligible')
            blocks = validate_maf(maf_path)
            files = [copy.deepcopy(training)] + [dict(kind='auxiliary',slot=None,path=relative,
                     sha256=manifest.sha256(source/relative)) for relative in (train+'.complete.json',chosen)]
        except (ValueError, OSError, StopIteration) as error:
            record['skipped'].append(f'{pair}: {error}'); continue
        for entry in files:
            _copy(target,source,entry['path'],entry['sha256'])
            existing[entry['path']] = entry
        record['pairs'][pair] = dict(train=train,maf=chosen,source_maf=chosen,maf_blocks=blocks,
                                    files=files,evidence=evidence,generated_by='lastal 1648')
    new['artifacts'] = list(existing.values())
    new['legacy_reuse'] = record
    manifest.atomic_json(target/manifest.MANIFEST,new)
    return record


def refresh_sorted(target):
    """Record maf-sort derivatives only after the caller's successful mafjoin stage."""
    value = manifest.load(target)
    changed = False
    for pair, record in value.get('legacy_reuse',{}).get('pairs',{}).items():
        original = manifest.artifact_path(target,record['maf'])
        if original.exists():
            if not pair_ready(target,pair): raise ValueError(f'Imported pair checksum mismatch: {pair}')
            continue
        if record['maf'].endswith('_sorted.maf'):
            raise ValueError(f'Imported sorted MAF missing: {pair}')
        sorted_name = record['maf'][:-4]+'_sorted.maf'
        sorted_path = manifest.artifact_path(target,sorted_name)
        validate_maf(sorted_path)
        for entry in record['files']:
            if entry['path'] != record['maf']: _verified_entry(target,entry)
        old_name = record['maf']
        item = dict(kind='auxiliary',slot=None,path=sorted_name,sha256=manifest.sha256(sorted_path))
        record['files'] = [a for a in record['files'] if a['path'] != old_name] + [item]
        record['maf'] = sorted_name
        record['transformation'] = dict(tool='maf-sort',source_path=old_name,
                                       source_sha256=next(a['sha256'] for a in value['artifacts'] if a['path']==old_name))
        value['artifacts'] = [a for a in value['artifacts'] if a['path'] not in (old_name,sorted_name)] + [item]
        changed = True
    if changed: manifest.atomic_json(Path(target)/manifest.MANIFEST,value)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    command=commands.add_parser('import')
    command.add_argument('target'); command.add_argument('source'); command.add_argument('log')
    for name in ('pair-ready','current-pair-ready'):
        command=commands.add_parser(name); command.add_argument('target'); command.add_argument('pair',choices=PAIRS)
    command=commands.add_parser('refresh-sorted'); command.add_argument('target')
    args=parser.parse_args()
    try:
        if args.command=='import':
            record=import_run(args.target,args.source,args.log)
            print('Certified legacy LAST 1648 pairs: '+(', '.join(record['pairs']) or 'none'))
            for reason in record['skipped']: print('Reuse skipped: '+reason,file=sys.stderr)
        elif args.command in ('pair-ready','current-pair-ready'):
            ready = pair_ready if args.command=='pair-ready' else current_pair_ready
            return 0 if ready(args.target,args.pair) else 1
        else: refresh_sorted(args.target)
    except (ValueError,OSError,KeyError) as error:
        print(f'Reuse error: {error}',file=sys.stderr); return 1
    return 0


if __name__=='__main__': sys.exit(main())
