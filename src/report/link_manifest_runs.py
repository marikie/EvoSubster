#!/usr/bin/env python3
import argparse
import os
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'lib'))
import run_manifest as runs


def main():
    parser = argparse.ArgumentParser(description='Link latest completed v2 runs by recorded trio identity')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('input_dirs', nargs='+')
    args = parser.parse_args()
    output = Path(args.output_dir).resolve()
    planned = {}
    for value in args.input_dirs:
        root = Path(value).resolve()
        if output == root:
            raise ValueError('Output cannot also be an input')
        for run in runs.discover(root):
            manifest = runs.validated(run, complete=True)
            key = manifest['trio_key']
            if key in planned and planned[key] != run:
                raise ValueError(f'Conflicting runs for trio: {key}')
            planned[key] = run
    for key, source in planned.items():
        target = output / key
        if target.exists() or target.is_symlink():
            if not target.is_symlink() or target.resolve() != source.resolve():
                raise ValueError(f'Destination collision: {target}')
    output.mkdir(parents=True, exist_ok=True)
    for key, source in planned.items():
        target = output / key
        if not target.is_symlink():
            target.symlink_to(os.path.relpath(source, output), target_is_directory=True)
    print(f'Trio candidates: {len(planned)}')


if __name__ == '__main__':
    main()
