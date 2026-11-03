#!/usr/bin/env python3
"""Parse LAST 1648 options and record settings before reusing pipeline output."""

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shlex
import sys


TOOLS = ('lastdb', 'lastal', 'last-train', 'last-split')
# Options taking a value, and switches. Positional operands belong to the pipeline.
SPECS = {
    'lastdb': ('uSRUwWQsmdai bBCP'.replace(' ', ''), 'pcqxDv', {'bits': True, 'circular': False}),
    'lastal': ('fDHErqpXabABcFxy zdelLmkWPKCsSiTnNRUu wGtgjJQ'.replace(' ', ''), 'v2M',
               {**{('split-' + x): True for x in 'fdctMSmsb'},
                'split': False, 'splice': False, 'split-n': False, 'reverse': False}),
    'last-train': ('rqpabABFDEsSC TRmkPXQ'.replace(' ', ''), 'v',
                   {**{x: True for x in ('fixmat', 'pid', 'postmask', 'sample-number', 'sample-length', 'scale')},
                    **{x: False for x in ('revsym', 'matsym', 'gapsym', 'codon', 'verbose')}}),
    'last-split': ('fgdctMSmsb', 'rnv',
                   {**{x: True for x in ('format', 'genome', 'direction', 'cis', 'trans', 'mean', 'sdev', 'mismap', 'score', 'bytes')},
                    **{x: False for x in ('reverse', 'no-split', 'verbose')}}),
}
PROTECTED = {
    'lastdb': {'-P', '-p', '-q', '-a', '-Q', '-x', '-D'},
    'lastal': {'-P', '-p', '-f', '-j', '-Q', '-2', '-F', '--split-f', '--split-n', '--splice', '--reverse'},
    'last-train': {'-P', '-Q', '-F', '--codon'},
    'last-split': {'-r', '--reverse', '-f', '--format', '-g', '--genome', '-n', '--no-split'},
}


def parse_options(tool, text):
    tokens = shlex.split(text)
    groups = []
    values, flags, long = SPECS[tool]
    index = 0
    while index < len(tokens):
        token = tokens[index]
        index += 1
        if token in ('-h', '--help', '-V', '--version', '--') or not token.startswith('-') or token == '-':
            raise ValueError('{}: operands, help and version flags are not allowed: {}'.format(tool, token))
        if token.startswith('--'):
            name, equal, value = token.partition('=')
            key = name[2:]
            if key not in long:
                raise ValueError('{}: unsupported LAST 1648 option: {}'.format(tool, name))
            if name in PROTECTED[tool]:
                raise ValueError('{}: {} is managed by the pipeline (use --thread for threads)'.format(tool, name))
            if long[key]:
                if not equal:
                    if index == len(tokens):
                        raise ValueError('{}: {} requires a value'.format(tool, name))
                    value = tokens[index]
                    index += 1
                if not value:
                    raise ValueError('{}: {} requires a non-empty value'.format(tool, name))
                groups.append([name, value])
            elif equal:
                raise ValueError('{}: {} does not accept a value'.format(tool, name))
            else:
                groups.append([name])
        else:
            chars = token[1:]
            while chars:
                key, chars = chars[0], chars[1:]
                name = '-' + key
                if name in PROTECTED[tool]:
                    raise ValueError('{}: {} is managed by the pipeline (use --thread for threads)'.format(tool, name))
                if key in values:
                    if chars:
                        value, chars = chars, ''
                    else:
                        if index == len(tokens):
                            raise ValueError('{}: {} requires a value'.format(tool, name))
                        value = tokens[index]
                        index += 1
                    if not value:
                        raise ValueError('{}: {} requires a non-empty value'.format(tool, name))
                    groups.append([name, value])
                elif key in flags:
                    groups.append([name])
                else:
                    raise ValueError('{}: unsupported LAST 1648 option: {}'.format(tool, name))
    return groups


def extras():
    data = json.loads(os.environ.get('EVOSUBSTER_LAST_OPTIONS', '{}'))
    if not isinstance(data, dict) or any(key not in TOOLS for key in data):
        raise ValueError('Invalid EVOSUBSTER_LAST_OPTIONS object')
    # The environment contains strings, never shell code or prevalidated argv.
    return {tool: parse_options(tool, data.get(tool, '')) for tool in TOOLS}


def merge(defaults, extra):
    overridden = {group[0] for group in extra}
    return [word for group in defaults if group[0] not in overridden for word in group] + [word for group in extra for word in group]


def config():
    custom = extras()
    threads = os.environ.get('THREAD_NUM', '8')
    if not threads.isdigit() or int(threads) < 1:
        raise ValueError('THREAD_NUM must be a positive integer')
    seed = os.environ.get('LASTDB_ALIGNMENT_SEED', 'RY128')
    custom_seed = any(group[0] == '-u' for group in custom['lastdb'])
    if not custom_seed and (not seed or any(not (c.isascii() and (c.isalnum() or c in '._-')) for c in seed)):
        raise ValueError('LASTDB_ALIGNMENT_SEED must not be empty or contain unsupported characters')
    return {
        'schema': 1,
        'lastdb_train': merge([['-P', threads], ['-c'], ['-u', 'RY4']], custom['lastdb']),
        'lastdb_align': merge([['-P', threads], ['-c'], ['-u', seed]], custom['lastdb']),
        'last_train': merge([['-P', threads], ['--revsym'], ['-C', '2']], custom['last-train']),
        'lastal': merge([['-P', threads], ['-i', os.environ.get('LASTAL_QUERY_BATCH_SIZE', '1M') or '1M'],
                         ['-H', '1'], ['-C', '2'], ['--split-f', 'MAF+'],
                         ['--split-b', os.environ.get('LASTAL_SPLIT_MEMORY', '8T') or '8T']], custom['lastal']),
        'last_split': ['-r'] + [word for group in custom['last-split'] for word in group],
    }


def database_label():
    words = config()['lastdb_align']
    seed = words[max(i for i, word in enumerate(words) if word == '-u') + 1]
    if seed and all(c.isascii() and (c.isalnum() or c in '._-') for c in seed):
        return seed
    # LAST also accepts seed-definition files. Never embed their paths in DB names.
    return 'custom-' + hashlib.sha256(json.dumps(words).encode()).hexdigest()[:16]


def guard(folder, stage):
    settings = config()
    if stage == 'training':
        settings = {key: settings[key] for key in ('schema', 'lastdb_train', 'last_train')}
        customized = bool(extras()['lastdb'] or extras()['last-train'])
    elif stage == 'database':
        settings = {key: settings[key] for key in ('schema', 'lastdb_align')}
        customized = bool(extras()['lastdb'])
    else:
        customized = any(extras().values())
    filename = {'run': 'metadata/last_options.json', 'training': 'last_training_options.json',
                'database': 'last_db_options.json', 'selection': 'last_options.json'}
    path = Path(str(folder) + '.last_options.json') if stage == 'alignment' else folder / filename[stage]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        guard_record(folder, stage, path, settings, customized)


def guard_record(folder, stage, path, settings, customized):
    if path.exists():
        if json.loads(path.read_text()) != settings:
            raise ValueError('LAST settings differ from existing results in {}. Use another date/--out-dir or --train-cache-dir.'.format(folder))
        return
    if customized and stage == 'alignment':
        if folder.exists() or folder.with_name(folder.stem + '_sorted.maf').exists():
            raise ValueError('Existing alignment has no LAST settings record: {}. Use another date/--out-dir.'.format(folder))
    elif customized and folder.exists():
        # Download metadata alone is not a completed analysis or reusable cache.
        existing = any(p.is_file() or p.is_symlink() for p in folder.rglob('*')
                       if p not in (path, path.with_suffix('.lock'))
                       and not (stage == 'run' and 'metadata' in p.relative_to(folder).parts))
        if existing:
            raise ValueError('Existing results/cache in {} have no LAST settings record; custom options cannot reuse them. Use another date/--out-dir or --train-cache-dir.'.format(folder))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix('.tmp.{}'.format(os.getpid()))
    temporary.write_text(json.dumps(settings, indent=2) + '\n')
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('validate', 'argv', 'guard', 'database-label'))
    parser.add_argument('value', nargs='?')
    parser.add_argument('stage', nargs='?', choices=('run', 'training', 'selection', 'alignment', 'database'))
    args = parser.parse_args()
    try:
        if args.action == 'validate':
            print(json.dumps(config(), sort_keys=True))
        elif args.action == 'database-label':
            print(database_label())
        elif args.action == 'argv':
            words = iter(config()[args.value])
            tool = {'lastdb_train': 'lastdb', 'lastdb_align': 'lastdb',
                    'last_train': 'last-train', 'last_split': 'last-split', 'lastal': 'lastal'}[args.value]
            for word in words:
                if word.startswith('--') and SPECS[tool][2].get(word[2:]):
                    word += '=' + next(words)
                elif len(word) == 2 and word[1] in SPECS[tool][0]:
                    value = next(words)
                    if word in ('-P', '-u'):
                        word += value
                    else:
                        sys.stdout.buffer.write(word.encode() + b'\0')
                        word = value
                sys.stdout.buffer.write(word.encode() + b'\0')
        else:
            guard(Path(args.value), args.stage)
    except (ValueError, KeyError, TypeError, OSError) as error:
        print('LAST options: ' + str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
