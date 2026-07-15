#!/usr/bin/env python3
"""Parse GFF features into deterministic disjoint strand-state segments."""

from __future__ import annotations

import gzip
import os
from bisect import bisect_right
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, TextIO, Tuple, Union


PathLike = Union[str, os.PathLike]


@dataclass(frozen=True)
class CDSSegment:
    """A zero-based, half-open genomic interval with one strand state."""

    start0: int
    end0: int
    strand: str

    def __post_init__(self) -> None:
        if not isinstance(self.start0, int) or not isinstance(self.end0, int):
            raise TypeError("segment coordinates must be integers")
        if self.start0 < 0 or self.end0 <= self.start0:
            raise ValueError("segment coordinates must satisfy 0 <= start0 < end0")
        if self.strand not in ("+", "-", "."):
            raise ValueError("segment strand must be '+', '-', or '.'")


class CDSIndex:
    """Per-seqid disjoint segments supporting binary-search point queries."""

    def __init__(self, segments_by_seqid: Mapping[str, Sequence[CDSSegment]]):
        ordered_seqids = tuple(sorted(segments_by_seqid))
        self._seqids = ordered_seqids
        self._segments_by_seqid = {
            seqid: tuple(segments_by_seqid[seqid]) for seqid in ordered_seqids
        }
        self._starts_by_seqid = {
            seqid: tuple(segment.start0 for segment in segments)
            for seqid, segments in self._segments_by_seqid.items()
        }

    @classmethod
    def from_gff(
        cls, path: PathLike, feature_type: str = "CDS"
    ) -> "CDSIndex":
        """Build an index from selected GFF3 features.

        GFF coordinates are converted from one-based inclusive coordinates to
        zero-based half-open coordinates. Duplicate and overlapping records are
        reduced to strand-state coverage rather than counted independently.
        """
        if not isinstance(feature_type, str) or not feature_type:
            raise ValueError("feature_type must be a non-empty string")

        events_by_seqid: Dict[str, Dict[int, List[int]]] = {}
        path_string = os.fspath(path)
        with _open_text(path_string) as stream:
            for line_number, raw_line in enumerate(stream, start=1):
                line = raw_line.rstrip("\r\n")
                if not line or line.startswith("#"):
                    continue

                fields = line.split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        "{} line {}: expected 9 tab-separated GFF fields".format(
                            path_string, line_number
                        )
                    )
                if fields[2] != feature_type:
                    continue

                seqid = fields[0]
                if not seqid or seqid == ".":
                    raise ValueError(
                        "{} line {}: selected feature has no seqid".format(
                            path_string, line_number
                        )
                    )
                start1, end1 = _parse_coordinates(
                    fields[3], fields[4], path_string, line_number
                )
                strand = fields[6]
                if strand not in ("+", "-"):
                    raise ValueError(
                        "{} line {}: selected feature strand must be '+' or '-'".format(
                            path_string, line_number
                        )
                    )

                start0 = start1 - 1
                end0 = end1
                strand_index = 0 if strand == "+" else 1
                seqid_events = events_by_seqid.setdefault(seqid, {})
                seqid_events.setdefault(start0, [0, 0])[strand_index] += 1
                seqid_events.setdefault(end0, [0, 0])[strand_index] -= 1

        segments_by_seqid = {
            seqid: _segments_from_events(events)
            for seqid, events in events_by_seqid.items()
        }
        return cls(segments_by_seqid)

    @property
    def seqids(self) -> Tuple[str, ...]:
        """Return indexed sequence identifiers in deterministic order."""
        return self._seqids

    def segments_for(self, seqid: str) -> Tuple[CDSSegment, ...]:
        """Return sorted disjoint segments for seqid, or an empty tuple."""
        return self._segments_by_seqid.get(seqid, ())

    def query(self, seqid: str, coord: int) -> Optional[CDSSegment]:
        """Return the segment containing a zero-based coordinate, if present."""
        if not isinstance(coord, int):
            raise TypeError("query coordinate must be an integer")
        if coord < 0:
            raise ValueError("query coordinate must be non-negative")

        segments = self._segments_by_seqid.get(seqid)
        if not segments:
            return None
        starts = self._starts_by_seqid[seqid]
        index = bisect_right(starts, coord) - 1
        if index >= 0 and coord < segments[index].end0:
            return segments[index]
        return None


def _open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def _parse_coordinates(
    start_text: str, end_text: str, path: str, line_number: int
) -> Tuple[int, int]:
    try:
        start1 = int(start_text)
        end1 = int(end_text)
    except ValueError as exc:
        raise ValueError(
            "{} line {}: selected feature coordinates must be integers".format(
                path, line_number
            )
        ) from exc
    if start1 < 1 or end1 < start1:
        raise ValueError(
            "{} line {}: selected feature coordinates must satisfy 1 <= start <= end".format(
                path, line_number
            )
        )
    return start1, end1


def _segments_from_events(events: Mapping[int, Sequence[int]]) -> Tuple[CDSSegment, ...]:
    segments: List[CDSSegment] = []
    plus_count = 0
    minus_count = 0
    previous_coord: Optional[int] = None

    for coord in sorted(events):
        if previous_coord is not None and previous_coord < coord:
            strand = _strand_state(plus_count, minus_count)
            if strand is not None:
                _append_coalesced(
                    segments, CDSSegment(previous_coord, coord, strand)
                )

        plus_delta, minus_delta = events[coord]
        plus_count += plus_delta
        minus_count += minus_delta
        if plus_count < 0 or minus_count < 0:
            raise ValueError("invalid interval event ordering")
        previous_coord = coord

    if plus_count != 0 or minus_count != 0:
        raise ValueError("unbalanced interval events")
    return tuple(segments)


def _strand_state(plus_count: int, minus_count: int) -> Optional[str]:
    if plus_count and minus_count:
        return "."
    if plus_count:
        return "+"
    if minus_count:
        return "-"
    return None


def _append_coalesced(
    segments: List[CDSSegment], segment: CDSSegment
) -> None:
    if (
        segments
        and segments[-1].end0 == segment.start0
        and segments[-1].strand == segment.strand
    ):
        previous = segments[-1]
        segments[-1] = CDSSegment(previous.start0, segment.end0, segment.strand)
        return
    segments.append(segment)


__all__ = ["CDSSegment", "CDSIndex"]
