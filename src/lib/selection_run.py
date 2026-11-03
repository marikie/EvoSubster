#!/usr/bin/env python3
"""Validate trio input and persist/resolve the settings needed to run it."""

import argparse
import csv
import json
import os
from pathlib import Path
import re
import sys
import tempfile

import last_options

ENV_DEFAULTS = {
    'LASTDB_ALIGNMENT_SEED': 'RY128',
    'LASTAL_QUERY_BATCH_SIZE': '1M',
    'LASTAL_SPLIT_MEMORY': '8T',
}
PATH_KEYS = ('genome_dir', 'out_dir', 'train_cache_dir')
ACC_COLUMNS = ('out_acc', 'in1_acc', 'in2_acc')


def read_trios(path):
    with open(path, newline='', encoding='utf-8-sig') as stream:
        reader = csv.reader(stream, delimiter='\t')
        header = next((row for row in reader if any(x.strip() for x in row)), None)
        # The selector historically writes a blank file when nothing is selected.
        if header is None:
            return []
        for col in ACC_COLUMNS:
            if header.count(col) != 1:
                raise ValueError('{}: expected exactly one {} column'.format(path, col))
        indices = [header.index(col) for col in ACC_COLUMNS]
        trios = []
        for row in reader:
            if not any(x.strip() for x in row):
                continue
            prefix = '{}: line {}'.format(path, reader.line_num)
            if len(row) != len(header):
                raise ValueError(prefix + ': column count differs from header')
            trio = [row[i].strip() for i in indices]
            if any(not re.fullmatch(r'GC[AF]_[0-9]+\.[0-9]+', acc) for acc in trio):
                raise ValueError(prefix + ': expected versioned GCA/GCF accessions')
            if len(set(trio)) != 3:
                raise ValueError(prefix + ': duplicate accession within trio')
            trios.append(trio)
        return trios


def validate_record(data):
    if not isinstance(data, dict) or data.get('schema') != 1:
        raise ValueError('unsupported selection_run schema (expected 1)')
    for key in ('date',) + PATH_KEYS:
        if not isinstance(data.get(key), str) or not data[key] or '\0' in data[key]:
            raise ValueError('selection_run: invalid ' + key)
    if any(not Path(data[key]).is_absolute() for key in PATH_KEYS):
        raise ValueError('selection_run paths must be absolute')
    if type(data.get('threads')) is not int or data['threads'] < 1:
        raise ValueError('selection_run: threads must be a positive integer')
    options = data.get('last_options')
    if not isinstance(options, dict) or set(options) != set(last_options.TOOLS):
        raise ValueError('selection_run: invalid last_options')
    for tool, value in options.items():
        if not isinstance(value, str) or '\0' in value:
            raise ValueError('selection_run: invalid LAST option string')
        last_options.parse_options(tool, value)
    env = data.get('last_environment')
    if not isinstance(env, dict) or set(env) != set(ENV_DEFAULTS):
        raise ValueError('selection_run: invalid last_environment')
    if any(not isinstance(value, str) or not value or '\0' in value for value in env.values()):
        raise ValueError('selection_run: invalid LAST environment value')
    return data


def save(args):
    options = json.loads(os.environ.get('EVOSUBSTER_LAST_OPTIONS', '{}'))
    data = dict(schema=1, date=args.date, threads=int(args.threads),
                last_options={tool: options.get(tool, '') for tool in last_options.TOOLS},
                last_environment={key: os.environ.get(key) or default
                                  for key, default in ENV_DEFAULTS.items()})
    data.update({key: str(Path(getattr(args, key)).resolve()) for key in PATH_KEYS})
    validate_record(data)
    os.environ['THREAD_NUM'] = str(data['threads'])
    last_options.config()
    path = Path(args.file)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, prefix='.selection-run-',
                                     delete=False) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(data, stream, indent=2)
            stream.write('\n')
            stream.close()
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)


def resolve(args):
    path = Path(args.config) if args.config else Path(args.trios).parent / 'selection_run.json'
    data = {}
    if args.config or path.exists() or path.is_symlink():
        data = validate_record(json.loads(path.read_text()))
    date = args.date or data.get('date', '')
    if not date:
        raise ValueError('DATE is required for a TSV without selection_run.json')
    paths = [getattr(args, key) or data.get(key, '') for key in PATH_KEYS]
    threads = args.threads or str(data.get('threads', 8))
    options = json.loads(os.environ.get('EVOSUBSTER_LAST_OPTIONS', '{}'))
    saved_options = data.get('last_options', {})
    options.update({key: value for key, value in saved_options.items()
                    if key not in args.explicit_last})
    env = {key: data.get('last_environment', {}).get(key, os.environ.get(key) or default)
           for key, default in ENV_DEFAULTS.items()}
    os.environ.update(env, THREAD_NUM=threads, EVOSUBSTER_LAST_OPTIONS=json.dumps(options))
    last_options.config()
    values = [date, *paths, threads, json.dumps(options), *env.values(),
              str(path.resolve()) if data else '(none; CLI/defaults)']
    for value in values:
        sys.stdout.buffer.write(value.encode() + b'\0')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    rows = commands.add_parser('trios')
    rows.add_argument('file')
    writer = commands.add_parser('save')
    writer.add_argument('file')
    resolver = commands.add_parser('resolve')
    resolver.add_argument('--trios', required=True)
    resolver.add_argument('--config', default='')
    resolver.add_argument('--explicit-last', action='append', default=[])
    for command in (writer, resolver):
        for key in ('date', 'threads') + PATH_KEYS:
            command.add_argument('--' + key.replace('_', '-'), default='',
                                 required=command is writer)
    args = parser.parse_args()
    try:
        if args.command == 'save':
            save(args)
        elif args.command == 'resolve':
            resolve(args)
        else:
            for trio in read_trios(args.file):
                print('\t'.join(trio))
    except (ValueError, TypeError, KeyError, OSError, csv.Error) as error:
        print('Trio input: ' + str(error), file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    sys.exit(main())
