#!/usr/bin/env python3

"""Remove terminal NCBI assembly accessions from Newick leaf labels."""

import argparse
import os
import re
import sys
import tempfile
from pathlib import Path
from typing import Tuple


ACCESSION_SUFFIX_RE = re.compile(r"_(?:GCA|GCF)_[0-9]+\.[0-9]+$")
LEAF_SEGMENT_END = set("():,;")


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


def rewrite_leaf_segment(segment: str) -> Tuple[str, bool, bool]:
    """Strip a logical leaf-label suffix while retaining comments and quoting."""
    logical_characters = []
    source_positions = []
    index = 0

    while index < len(segment):
        char = segment[index]
        if char == "[":
            index = consume_comment(segment, index)
            continue
        if char == "'":
            end = consume_quoted_label(segment, index)
            quoted_index = index + 1
            while quoted_index < end - 1:
                logical_characters.append(segment[quoted_index])
                source_positions.append(quoted_index)
                if (
                    segment[quoted_index] == "'"
                    and quoted_index + 1 < end - 1
                    and segment[quoted_index + 1] == "'"
                ):
                    quoted_index += 2
                else:
                    quoted_index += 1
            index = end
            continue
        if not char.isspace():
            logical_characters.append(char)
            source_positions.append(index)
        index += 1

    logical_label = "".join(logical_characters)
    match = ACCESSION_SUFFIX_RE.search(logical_label)
    if match is None:
        return segment, False, bool(logical_label)

    removed_positions = set(source_positions[match.start() : match.end()])
    rewritten = "".join(
        char for position, char in enumerate(segment) if position not in removed_positions
    )
    return rewritten, True, True


def consume_leaf_segment(text: str, start: int) -> int:
    """Return the end of a leaf-label region, including embedded comments."""
    index = start
    while index < len(text):
        char = text[index]
        if char in LEAF_SEGMENT_END:
            break
        if char == "[":
            index = consume_comment(text, index)
            continue
        if char == "'":
            index = consume_quoted_label(text, index)
            continue
        index += 1
    return index


def rewrite_newick(text: str) -> Tuple[str, int, int]:
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
            end = consume_leaf_segment(text, index)
            segment, changed, has_label = rewrite_leaf_segment(text[index:end])
            output.append(segment)
            if has_label:
                leaf_count += 1
                converted += int(changed)
                expecting_leaf = False
            if end > index:
                index = end
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

    return "".join(output), converted, leaf_count


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

    if input_path.resolve(strict=False) == output_path.resolve(strict=False):
        print("error: --input and --output must be different paths.", file=sys.stderr)
        return 2
    if not input_path.is_file():
        print(f"error: input Newick file does not exist: {input_path}", file=sys.stderr)
        return 2
    if not output_path.parent.is_dir():
        print(f"error: output directory does not exist: {output_path.parent}", file=sys.stderr)
        return 2
    try:
        rewritten, converted, leaf_count = rewrite_newick(read_text(input_path))
        write_text_atomic(output_path, rewritten, overwrite=args.force)
    except FileExistsError:
        print(
            f"error: output file already exists: {output_path} (use --force to overwrite)",
            file=sys.stderr,
        )
        return 2
    except (OSError, UnicodeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    unchanged = leaf_count - converted
    print(
        f"Converted {converted} of {leaf_count} leaf labels; output: {output_path}",
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
