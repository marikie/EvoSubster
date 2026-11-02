#!/usr/bin/env python3
"""Calculate CDS union length / whole FASTA length (including N) * 100."""

import argparse
from collections import defaultdict
import os
from pathlib import Path
import sys
import tempfile


def fasta_lengths(path):
    """Stream sequences, retaining only lengths keyed by the first header token."""
    lengths = {}
    seq_id = None
    with open(path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields or fields[0] in lengths:
                    raise ValueError(f"{path}:{line_number}: empty or duplicate FASTA ID")
                seq_id = fields[0]
                lengths[seq_id] = 0
            elif seq_id is None:
                raise ValueError(f"{path}:{line_number}: sequence before FASTA header")
            else:
                lengths[seq_id] += len("".join(line.split()))
    if not lengths or any(length == 0 for length in lengths.values()):
        raise ValueError(f"{path}: empty FASTA or empty sequence")
    return lengths


def coding_percentage(fasta_path, gff_path):
    lengths = fasta_lengths(fasta_path)
    intervals = defaultdict(list)
    with open(gff_path, encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.strip() == "##FASTA":
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            location = f"{gff_path}:{line_number}"
            if len(fields) != 9:
                raise ValueError(f"{location}: expected 9 GFF columns")
            if fields[2] != "CDS":
                continue
            seq_id = fields[0]
            if seq_id not in lengths:
                raise ValueError(f"{location}: CDS sequence ID {seq_id!r} absent from FASTA")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"{location}: CDS coordinates must be integers") from exc
            if not 1 <= start <= end <= lengths[seq_id]:
                raise ValueError(f"{location}: CDS coordinates outside FASTA sequence bounds")
            intervals[seq_id].append((start - 1, end))

    coding_bases = 0
    for ranges in intervals.values():
        merged_start, merged_end = 0, 0
        for start, end in sorted(ranges):
            if start > merged_end:
                coding_bases += merged_end - merged_start
                merged_start, merged_end = start, end
            else:
                merged_end = max(merged_end, end)
        coding_bases += merged_end - merged_start
    return 100 * coding_bases / sum(lengths.values())


def write_atomic(path, content):
    """Expose a completed result only after the entire write succeeds."""
    path = Path(path)
    staged = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        ) as handle:
            staged = Path(handle.name)
            handle.write(content)
        os.replace(staged, path)
    finally:
        if staged is not None and staged.exists():
            staged.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fasta", help="Whole genome FASTA (including unannotated contigs)")
    parser.add_argument("gff", help="GFF annotation matching the FASTA sequence IDs")
    parser.add_argument("--output", type=Path, help="Atomically write result to this .out file")
    args = parser.parse_args()
    try:
        percentage = coding_percentage(args.fasta, args.gff)
        result = f"Coding region percentage: {percentage:.2f}%\n"
        if args.output is None:
            sys.stdout.write(result)
        else:
            write_atomic(args.output, result)
    except (OSError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
