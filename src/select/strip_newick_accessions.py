#!/usr/bin/env python3

"""Remove terminal NCBI assembly accessions from Newick leaf labels."""

import argparse
import os
import re
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple


ACCESSION_SUFFIX_RE = re.compile(r"_(?:GCA|GCF)_[0-9]+\.[0-9]+$")
LEAF_SEGMENT_END = set("():,;")


@dataclass(frozen=True)
class LeafScan:
    end: int
    logical_label: str
    source_positions: Tuple[int, ...]
    has_label: bool


@dataclass(frozen=True)
class ConversionResult:
    text: str
    converted: int
    leaf_count: int


class CliUsageError(ValueError):
    """Invalid CLI path configuration that retains exit code 2."""


def consume_comment(text: str, start: int) -> int:
    """Return the first index after a bracketed Newick comment."""
    end = text.find("]", start + 1)
    if end < 0:
        raise ValueError("unterminated Newick comment")
    return end + 1


def consume_quoted_label(text: str, start: int) -> int:
    """Return the first index after a single-quoted Newick label."""
    index = start + 1
    while index < len(text):
        if text[index] != "'":
            index += 1
            continue
        if index + 1 < len(text) and text[index + 1] == "'":
            index += 2
            continue
        return index + 1
    raise ValueError("unterminated quoted Newick label")


def scan_leaf_segment(text: str, start: int) -> LeafScan:
    """Scan a leaf-label region and retain its logical source positions."""
    logical_characters = []
    source_positions = []
    index = start

    while index < len(text):
        char = text[index]
        if char in LEAF_SEGMENT_END:
            break
        if char == "[":
            index += 1
            while index < len(text) and text[index] != "]":
                index += 1
            if index == len(text):
                raise ValueError("unterminated Newick comment")
            index += 1
            continue
        if char == "'":
            index += 1
            while index < len(text):
                quoted_char = text[index]
                if quoted_char != "'":
                    logical_characters.append(quoted_char)
                    source_positions.append(index)
                    index += 1
                    continue
                if (
                    index + 1 < len(text)
                    and text[index + 1] == "'"
                ):
                    logical_characters.append(quoted_char)
                    source_positions.append(index)
                    index += 2
                    continue
                index += 1
                break
            else:
                raise ValueError("unterminated quoted Newick label")
            continue
        if not char.isspace():
            logical_characters.append(char)
            source_positions.append(index)
        index += 1

    logical_label = "".join(logical_characters)
    return LeafScan(index, logical_label, tuple(source_positions), bool(logical_label))


def rewrite_leaf_segment(text: str, start: int, scan: LeafScan) -> Tuple[str, bool]:
    """Strip a scanned logical leaf-label suffix while retaining its source text."""
    match = ACCESSION_SUFFIX_RE.search(scan.logical_label)
    if match is None:
        return text[start : scan.end], False

    removed_positions = set(scan.source_positions[match.start() : match.end()])
    rewritten = "".join(
        char
        for position, char in enumerate(text[start : scan.end], start)
        if position not in removed_positions
    )
    return rewritten, True


def rewrite_newick(text: str) -> ConversionResult:
    """Rewrite leaf labels while preserving all other Newick text verbatim."""
    output = []
    index = 0
    expecting_leaf = True
    converted = 0
    leaf_count = 0

    while index < len(text):
        char = text[index]

        if char == "(":
            output.append(char)
            expecting_leaf = True
            index += 1
            continue

        if char == ",":
            output.append(char)
            expecting_leaf = True
            index += 1
            continue

        if char == ")":
            output.append(char)
            expecting_leaf = False
            index += 1
            continue

        if expecting_leaf:
            scan = scan_leaf_segment(text, index)
            segment, changed = rewrite_leaf_segment(text, index, scan)
            output.append(segment)
            if scan.has_label:
                leaf_count += 1
                converted += int(changed)
                expecting_leaf = False
            if scan.end > index:
                index = scan.end
                continue

        if char == "[":
            end = consume_comment(text, index)
            output.append(text[index:end])
            index = end
            continue

        if char == "'":
            end = consume_quoted_label(text, index)
            output.append(text[index:end])
            index = end
            continue

        output.append(char)
        index += 1

    return ConversionResult("".join(output), converted, leaf_count)


def read_text(path: Path) -> str:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return handle.read()


def write_text_atomic(path: Path, text: str, overwrite: bool) -> None:
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            "w",
            encoding="utf-8",
            newline="",
            dir=path.parent,
            prefix=f".{path.name}.",
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            handle.write(text)
        if overwrite:
            os.replace(temporary_path, path)
        else:
            os.link(temporary_path, path)
            temporary_path.unlink()
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def convert_file(input_path: Path, output_path: Path, overwrite: bool) -> ConversionResult:
    """Validate, convert, and atomically publish a Newick file."""
    if input_path.resolve(strict=False) == output_path.resolve(strict=False):
        raise CliUsageError("--input and --output must be different paths.")
    if not input_path.is_file():
        raise CliUsageError(f"input Newick file does not exist: {input_path}")
    if not output_path.parent.is_dir():
        raise CliUsageError(f"output directory does not exist: {output_path.parent}")

    result = rewrite_newick(read_text(input_path))
    try:
        write_text_atomic(output_path, result.text, overwrite=overwrite)
    except FileExistsError as exc:
        raise CliUsageError(
            f"output file already exists: {output_path} (use --force to overwrite)"
        ) from exc
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Remove a terminal _GCA_...version or _GCF_...version suffix from each "
            "Newick leaf label without changing the rest of the tree text."
        )
    )
    parser.add_argument("--input", required=True, type=Path, help="Input Newick file.")
    parser.add_argument("--output", required=True, type=Path, help="Output Newick file.")
    parser.add_argument(
        "--force", action="store_true", help="Overwrite an existing output file."
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_path = args.input.expanduser()
    output_path = args.output.expanduser()

    try:
        result = convert_file(input_path, output_path, overwrite=args.force)
    except CliUsageError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except (OSError, UnicodeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    unchanged = result.leaf_count - result.converted
    print(
        f"Converted {result.converted} of {result.leaf_count} leaf labels; output: {output_path}",
        file=sys.stderr,
    )
    if unchanged:
        noun = "label" if unchanged == 1 else "labels"
        print(
            f"warning: {unchanged} leaf {noun} had no terminal accession and was left unchanged.",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
