#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
from typing import Iterable

PATTERN = "*_ncds.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create symlinks under the output directory for all matching TSV files found recursively in the input directory."
    )
    parser.add_argument("input_dir", help="Root directory to search for TSV files.")
    parser.add_argument("output_dir", help="Directory where symlinks will be placed.")
    parser.add_argument(
        "--pattern",
        default=PATTERN,
        help=f"Glob pattern to match files (default: {PATTERN}).",
    )
    return parser.parse_args()


def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def iter_matching_files(root: Path, pattern: str) -> Iterable[Path]:
    return root.rglob(pattern)


def create_symlink(src: Path, dest: Path) -> None:
    dest.symlink_to(src)


def process_links(input_dir: Path, output_dir: Path, pattern: str) -> None:
    input_dir = input_dir.resolve()
    output_dir = output_dir.resolve()
    if input_dir == output_dir:
        raise ValueError("Input and output directories must be different")
    planned = {}
    for src in sorted(iter_matching_files(input_dir, pattern)):
        if output_dir == src or output_dir in src.parents:
            continue
        if not src.is_file():
            continue
        source = src.resolve()
        dest = output_dir / src.name
        if dest in planned and planned[dest] != source:
            raise ValueError(f"Conflicting source basename: {src.name}")
        if dest.exists() or dest.is_symlink():
            if not dest.is_symlink() or dest.resolve() != source:
                raise ValueError(f"Destination already belongs to another file: {dest}")
        planned[dest] = source
    ensure_directory(output_dir)
    created = 0
    for dest, source in planned.items():
        if dest.is_symlink():
            continue
        create_symlink(source, dest)
        created += 1
    print(f"New symlinks created: {created}")


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir).resolve()
    output_dir = Path(args.output_dir).resolve()

    if not input_dir.is_dir():
        print(f"Error: input directory does not exist: {input_dir}", file=sys.stderr)
        sys.exit(1)

    process_links(input_dir, output_dir, args.pattern)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

