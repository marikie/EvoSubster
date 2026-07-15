#!/usr/bin/env python3
"""Strict, streaming MAF I/O and forward-coordinate helpers."""

from dataclasses import dataclass
import gzip
import os
import sys
from typing import Iterable, Iterator, List, Optional, TextIO, Union


Path = Union[str, os.PathLike]


@dataclass(frozen=True)
class MafRow:
    """One MAF ``s`` row with the source-relative start from the file."""

    seqid: str
    start0: int
    span: int
    strand: str
    seqlen: int
    sequence: str


def open_text(path: Path) -> TextIO:
    """Open a plain or gzip-compressed text file, or return stdin for ``-``."""
    path_string = os.fspath(path)
    if path_string == "-":
        return sys.stdin
    if path_string.endswith(".gz"):
        return gzip.open(path_string, "rt")
    return open(path_string)


def _line_error(line_number: int, message: str) -> ValueError:
    return ValueError("MAF line {}: {}".format(line_number, message))


def _parse_row(fields: List[str], line_number: int) -> MafRow:
    if len(fields) != 7:
        raise _line_error(
            line_number,
            "expected 7 fields in sequence row, found {}".format(len(fields)),
        )

    _marker, seqid, start_text, span_text, strand, seqlen_text, sequence = fields
    try:
        start0 = int(start_text)
        span = int(span_text)
        seqlen = int(seqlen_text)
    except ValueError as exc:
        raise _line_error(line_number, "start, span, and seqlen must be integers") from exc

    if not seqid:
        raise _line_error(line_number, "seqid must not be empty")
    if start0 < 0 or span < 0 or seqlen < 0:
        raise _line_error(line_number, "start, span, and seqlen must be non-negative")
    if strand not in ("+", "-"):
        raise _line_error(line_number, "strand must be '+' or '-'")
    if start0 + span > seqlen:
        raise _line_error(line_number, "row coordinates exceed the source sequence")

    observed_span = sum(base != "-" for base in sequence)
    if observed_span != span:
        raise _line_error(
            line_number,
            "declared span {} does not match {} non-gap bases".format(
                span, observed_span
            ),
        )

    return MafRow(seqid, start0, span, strand, seqlen, sequence)


def _validated_block(rows: List[MafRow], block_number: int) -> List[MafRow]:
    aligned_lengths = {len(row.sequence) for row in rows}
    if len(aligned_lengths) != 1:
        raise ValueError(
            "MAF block {}: sequence rows have unequal aligned lengths {}".format(
                block_number, sorted(aligned_lengths)
            )
        )
    return rows


def iter_maf_blocks(lines: Iterable[str]) -> Iterator[List[MafRow]]:
    """Yield validated sequence rows from each non-empty MAF block."""
    rows = []
    block_number = 0

    for line_number, line in enumerate(lines, 1):
        fields = line.split()
        if not fields:
            if rows:
                block_number += 1
                yield _validated_block(rows, block_number)
                rows = []
            continue

        marker = fields[0]
        if marker == "a":
            if rows:
                block_number += 1
                yield _validated_block(rows, block_number)
                rows = []
        elif marker == "s":
            rows.append(_parse_row(fields, line_number))

    if rows:
        block_number += 1
        yield _validated_block(rows, block_number)


def forward_coord(row: MafRow, ungapped_offset: int) -> int:
    """Convert an ungapped row offset to a forward-strand source coordinate."""
    if row.strand not in ("+", "-"):
        raise ValueError("MAF row strand must be '+' or '-'")
    if not 0 <= ungapped_offset < row.span:
        raise ValueError(
            "ungapped offset {} is outside row span {}".format(
                ungapped_offset, row.span
            )
        )

    source_offset = row.start0 + ungapped_offset
    if row.strand == "+":
        return source_offset
    return row.seqlen - 1 - source_offset


def aligned_forward_coords(row: MafRow) -> List[Optional[int]]:
    """Return one forward coordinate per aligned column, using None for gaps."""
    observed_span = sum(base != "-" for base in row.sequence)
    if observed_span != row.span:
        raise ValueError(
            "MAF row span {} does not match {} non-gap bases".format(
                row.span, observed_span
            )
        )

    coordinates = []
    ungapped_offset = 0
    for base in row.sequence:
        if base == "-":
            coordinates.append(None)
        else:
            coordinates.append(forward_coord(row, ungapped_offset))
            ungapped_offset += 1
    return coordinates
