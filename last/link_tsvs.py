#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
from typing import Iterable

PATTERN = "*_maflinked_ncds.tsv"


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
    ensure_directory(output_dir)

    total_found = 0
    created_links = 0
    replaced_links = 0
    skipped_entries = 0

    for src in iter_matching_files(input_dir, pattern):
        total_found += 1
        dest = output_dir / src.name

        if dest.exists() or dest.is_symlink():
            if dest.is_dir() and not dest.is_symlink():
                print(
                    f"Warning: destination is a directory, skipping: {dest}",
                    file=sys.stderr,
                )
                skipped_entries += 1
                continue
            try:
                dest.unlink()
                replaced_links += 1
            except OSError as exc:
                print(
                    f"Warning: failed to remove existing entry {dest}: {exc}",
                    file=sys.stderr,
                )
                skipped_entries += 1
                continue

        create_symlink(src, dest)
        created_links += 1

    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Pattern: {pattern}")
    print(f"Total matches found: {total_found}")
    print(f"New symlinks created: {created_links}")
    print(f"Existing entries replaced: {replaced_links}")
    print(f"Entries skipped: {skipped_entries}")


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

