#!/usr/bin/env python3
"""Independent CDS annotation and strand consistency QC."""

import argparse
import csv
import bisect
import gzip
import hashlib
import os
import sys
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple


ALIGN_DIR = Path(__file__).resolve().parents[1] / "align"
if str(ALIGN_DIR) not in sys.path:
    sys.path.insert(0, str(ALIGN_DIR))

from cds_intervals import CDSIndex  # noqa: E402
from maf_io import aligned_forward_coords, iter_maf_blocks, open_text  # noqa: E402


CLASS_NONCODING = "nonCodingIncluded"
CLASS_CONSISTENT = "codingStrandConsistent"
CLASS_CONTRADICTION = "strandContradictions"
CLASS_INSUFFICIENT = "insufficientEvidence"
CLASSES = (
    CLASS_NONCODING,
    CLASS_CONSISTENT,
    CLASS_CONTRADICTION,
    CLASS_INSUFFICIENT,
)

RAW_COLUMNS = [
    "trio",
    "block_index",
    "run_index",
    "run_col_start0",
    "run_col_end0",
    "a_seqid",
    "a_start0",
    "a_end0",
    "a_coord_first0",
    "a_coord_last0",
    "a_cds_strand",
    "a_maf_strand",
    "a_segment_id",
    "b_seqid",
    "b_coord_min",
    "b_coord_max",
    "b_coord_first0",
    "b_coord_last0",
    "b_non_gap_bases",
    "b_cds_covered_bases",
    "b_strand_conflict_bases",
    "b_cov_frac",
    "b_cds_strand",
    "b_maf_strand",
    "c_seqid",
    "c_coord_min",
    "c_coord_max",
    "c_coord_first0",
    "c_coord_last0",
    "c_non_gap_bases",
    "c_cds_covered_bases",
    "c_strand_conflict_bases",
    "c_cov_frac",
    "c_cds_strand",
    "c_maf_strand",
    "match_a",
    "match_b",
    "match_c",
    "classification",
    "classification_reason",
]

SUMMARY_COLUMNS = [
    "trio",
    "min_non_gap",
    CLASS_NONCODING,
    CLASS_CONSISTENT,
    CLASS_CONTRADICTION,
    CLASS_INSUFFICIENT,
    "total_runs",
    "classified_runs",
    "blocks_seen",
    "blocks_lt3",
]

DEFAULT_TRIOS = ("lenEdo", "puffer", "styPis")
MANIFEST_COLUMNS = (
    "trio",
    "role",
    "accession",
    "annotation_source",
    "package_filename",
    "bytes",
    "sha256",
    "url",
)


@dataclass(frozen=True)
class ClassificationResult:
    classification: str
    reason: str
    b_cds_strand: str = "."
    c_cds_strand: str = "."
    match_a: Optional[bool] = None
    match_b: Optional[bool] = None
    match_c: Optional[bool] = None


@dataclass(frozen=True)
class SpeciesEvidence:
    non_gap: int
    covered: int
    plus: int
    minus: int
    conflict: int
    coord_min: Optional[int]
    coord_max: Optional[int]
    coord_first: Optional[int]
    coord_last: Optional[int]


@dataclass(frozen=True)
class AlignmentRun:
    block_index: int
    run_index: int
    col_start0: int
    col_end0: int
    segment: Any
    a_coord_first: int
    a_coord_last: int


@dataclass(frozen=True)
class TrioAnalysis:
    summary: Dict[str, Any]
    sensitivity: Tuple[Dict[str, Any], ...]


def _majority_strand(plus: int, minus: int) -> str:
    if plus == minus:
        return "."
    return "+" if plus > minus else "-"


def _validate_counts(non_gap: int, covered: int, plus: int, minus: int,
                     conflict: int, label: str) -> None:
    values = (non_gap, covered, plus, minus, conflict)
    if any(value < 0 for value in values):
        raise ValueError("{} evidence counts must be non-negative".format(label))
    if covered > non_gap:
        raise ValueError("{} CDS-covered count exceeds non-gap count".format(label))
    if plus + minus + conflict > covered:
        raise ValueError("{} strand-state counts exceed CDS-covered count".format(label))


def classify_evidence(
    *,
    a_cds_strand: str,
    a_maf_strand: str,
    b_maf_strand: str,
    c_maf_strand: str,
    b_non_gap: int,
    b_covered: int,
    b_plus: int,
    b_minus: int,
    b_conflict: int,
    c_non_gap: int,
    c_covered: int,
    c_plus: int,
    c_minus: int,
    c_conflict: int,
    min_non_gap: int,
) -> ClassificationResult:
    """Classify one A-CDS alignment run from precomputed evidence counts."""
    if min_non_gap < 1:
        raise ValueError("min_non_gap must be at least 1")
    for value, label in (
        (a_cds_strand, "A CDS"),
        (a_maf_strand, "A MAF"),
        (b_maf_strand, "B MAF"),
        (c_maf_strand, "C MAF"),
    ):
        allowed = ("+", "-", ".") if label == "A CDS" else ("+", "-")
        if value not in allowed:
            raise ValueError("{} strand must be one of {}".format(label, allowed))

    _validate_counts(b_non_gap, b_covered, b_plus, b_minus, b_conflict, "B")
    _validate_counts(c_non_gap, c_covered, c_plus, c_minus, c_conflict, "C")

    if a_cds_strand == ".":
        return ClassificationResult(CLASS_INSUFFICIENT,
                                    "a_cds_strand_ambiguous")
    if b_non_gap < min_non_gap or c_non_gap < min_non_gap:
        return ClassificationResult(CLASS_INSUFFICIENT,
                                    "b_or_c_below_min_non_gap")

    b_cov = b_covered / b_non_gap
    c_cov = c_covered / c_non_gap
    if b_cov < 0.5 or c_cov < 0.5:
        return ClassificationResult(CLASS_NONCODING,
                                    "b_or_c_below_cds_coverage")

    if b_conflict or c_conflict:
        return ClassificationResult(CLASS_INSUFFICIENT,
                                    "b_or_c_cds_strand_conflict")

    b_strand = _majority_strand(b_plus, b_minus)
    c_strand = _majority_strand(c_plus, c_minus)
    if b_strand == "." or c_strand == ".":
        return ClassificationResult(
            CLASS_INSUFFICIENT,
            "b_or_c_cds_strand_tie",
            b_cds_strand=b_strand,
            c_cds_strand=c_strand,
        )

    match_a = a_cds_strand == a_maf_strand
    match_b = b_strand == b_maf_strand
    match_c = c_strand == c_maf_strand
    classification = (
        CLASS_CONSISTENT
        if match_a == match_b == match_c
        else CLASS_CONTRADICTION
    )
    reason = (
        "strand_matches_equal"
        if classification == CLASS_CONSISTENT
        else "strand_matches_mixed"
    )
    return ClassificationResult(
        classification,
        reason,
        b_cds_strand=b_strand,
        c_cds_strand=c_strand,
        match_a=match_a,
        match_b=match_b,
        match_c=match_c,
    )


def annotate_monotonic(index: Any, seqid: str,
                       coords: Sequence[Optional[int]]) -> List[Any]:
    """Map monotonic forward coordinates to disjoint CDS segments in O(n)."""
    segments = tuple(index.segments_for(seqid))
    result = [None] * len(coords)
    non_gap = [(i, coord) for i, coord in enumerate(coords) if coord is not None]
    if not segments or not non_gap:
        return result

    values = [coord for _i, coord in non_gap]
    if len(values) > 1:
        deltas = [right - left for left, right in zip(values, values[1:])]
        if any(delta == 0 for delta in deltas):
            raise ValueError("non-gap MAF coordinates must be strictly monotonic")
        increasing = deltas[0] > 0
        if any((delta > 0) != increasing for delta in deltas[1:]):
            raise ValueError("non-gap MAF coordinates change direction")
    else:
        increasing = True

    starts = [segment.start0 for segment in segments]
    first = values[0]
    pos = bisect.bisect_right(starts, first) - 1
    if increasing:
        pos = max(pos, 0)
        for out_index, coord in non_gap:
            while pos < len(segments) and segments[pos].end0 <= coord:
                pos += 1
            if pos >= len(segments):
                break
            segment = segments[pos]
            if segment.start0 <= coord < segment.end0:
                result[out_index] = segment
    else:
        for out_index, coord in non_gap:
            while pos >= 0 and segments[pos].start0 > coord:
                pos -= 1
            if pos < 0:
                break
            segment = segments[pos]
            if segment.start0 <= coord < segment.end0:
                result[out_index] = segment
    return result


def collect_species_evidence(index: Any, seqid: str,
                             coords: Sequence[Optional[int]]) -> SpeciesEvidence:
    """Summarize independent CDS coverage and strand states for one run."""
    segments = annotate_monotonic(index, seqid, coords)
    observed = [coord for coord in coords if coord is not None]
    covered = plus = minus = conflict = 0
    for coord, segment in zip(coords, segments):
        if coord is None or segment is None:
            continue
        covered += 1
        if segment.strand == "+":
            plus += 1
        elif segment.strand == "-":
            minus += 1
        elif segment.strand == ".":
            conflict += 1
        else:
            raise ValueError("unexpected CDS segment strand: {!r}".format(
                segment.strand))
    return SpeciesEvidence(
        non_gap=len(observed),
        covered=covered,
        plus=plus,
        minus=minus,
        conflict=conflict,
        coord_min=min(observed) if observed else None,
        coord_max=max(observed) if observed else None,
        coord_first=observed[0] if observed else None,
        coord_last=observed[-1] if observed else None,
    )


def _optional(value: Optional[int]) -> Any:
    return "" if value is None else value


def _format_match(value: Optional[bool]) -> str:
    if value is None:
        return ""
    return "true" if value else "false"


def _evidence_strand(evidence: SpeciesEvidence) -> str:
    if evidence.conflict:
        return "."
    return _majority_strand(evidence.plus, evidence.minus)


def _summary_row(trio: str, threshold: int, counts: Dict[str, int],
                 blocks_seen: int, blocks_lt3: int) -> Dict[str, Any]:
    total = sum(counts.get(label, 0) for label in CLASSES)
    return {
        "trio": trio,
        "min_non_gap": threshold,
        CLASS_NONCODING: counts.get(CLASS_NONCODING, 0),
        CLASS_CONSISTENT: counts.get(CLASS_CONSISTENT, 0),
        CLASS_CONTRADICTION: counts.get(CLASS_CONTRADICTION, 0),
        CLASS_INSUFFICIENT: counts.get(CLASS_INSUFFICIENT, 0),
        "total_runs": total,
        "classified_runs": total - counts.get(CLASS_INSUFFICIENT, 0),
        "blocks_seen": blocks_seen,
        "blocks_lt3": blocks_lt3,
    }


def analyze_trio(
    trio: str,
    maf_path: Path,
    a_gff_path: Path,
    b_gff_path: Path,
    c_gff_path: Path,
    raw_path: Path,
    min_non_gap: int = 10,
    sensitivity_thresholds: Sequence[int] = (1, 10, 30),
) -> TrioAnalysis:
    """Stream one trio and write one raw evidence row per A-CDS run."""
    thresholds = tuple(sorted(set(sensitivity_thresholds) | {min_non_gap}))
    if not thresholds or any(threshold < 1 for threshold in thresholds):
        raise ValueError("sensitivity thresholds must be positive integers")

    a_index = CDSIndex.from_gff(a_gff_path)
    b_index = CDSIndex.from_gff(b_gff_path)
    c_index = CDSIndex.from_gff(c_gff_path)
    counts = {threshold: {label: 0 for label in CLASSES}
              for threshold in thresholds}
    blocks_seen = 0
    blocks_lt3 = 0
    raw_rows_written = 0

    with atomic_csv(raw_path, RAW_COLUMNS) as raw_writer:
        with open_text(maf_path) as maf_handle:
            for block_index, rows in enumerate(iter_maf_blocks(maf_handle)):
                blocks_seen += 1
                if len(rows) < 3:
                    blocks_lt3 += 1
                    continue
                if len(rows) != 3:
                    raise ValueError(
                        "MAF block {} must contain exactly 3 sequence rows; found {}"
                        .format(block_index, len(rows))
                    )

                a_row, b_row, c_row = rows
                a_coords = aligned_forward_coords(a_row)
                b_coords = aligned_forward_coords(b_row)
                c_coords = aligned_forward_coords(c_row)
                a_segments = annotate_monotonic(a_index, a_row.seqid, a_coords)

                for run in iter_a_runs(block_index, a_coords, a_segments):
                    run_slice = slice(run.col_start0, run.col_end0)
                    b_evidence = collect_species_evidence(
                        b_index, b_row.seqid, b_coords[run_slice]
                    )
                    c_evidence = collect_species_evidence(
                        c_index, c_row.seqid, c_coords[run_slice]
                    )
                    kwargs = {
                        "a_cds_strand": run.segment.strand,
                        "a_maf_strand": a_row.strand,
                        "b_maf_strand": b_row.strand,
                        "c_maf_strand": c_row.strand,
                        "b_non_gap": b_evidence.non_gap,
                        "b_covered": b_evidence.covered,
                        "b_plus": b_evidence.plus,
                        "b_minus": b_evidence.minus,
                        "b_conflict": b_evidence.conflict,
                        "c_non_gap": c_evidence.non_gap,
                        "c_covered": c_evidence.covered,
                        "c_plus": c_evidence.plus,
                        "c_minus": c_evidence.minus,
                        "c_conflict": c_evidence.conflict,
                    }
                    by_threshold = {}
                    for threshold in thresholds:
                        classified = classify_evidence(
                            **kwargs, min_non_gap=threshold
                        )
                        by_threshold[threshold] = classified
                        counts[threshold][classified.classification] += 1

                    primary = by_threshold[min_non_gap]
                    b_cov = (b_evidence.covered / b_evidence.non_gap
                             if b_evidence.non_gap else 0.0)
                    c_cov = (c_evidence.covered / c_evidence.non_gap
                             if c_evidence.non_gap else 0.0)
                    raw_writer.writerow({
                        "trio": trio,
                        "block_index": block_index,
                        "run_index": run.run_index,
                        "run_col_start0": run.col_start0,
                        "run_col_end0": run.col_end0,
                        "a_seqid": a_row.seqid,
                        "a_start0": run.segment.start0,
                        "a_end0": run.segment.end0,
                        "a_coord_first0": run.a_coord_first,
                        "a_coord_last0": run.a_coord_last,
                        "a_cds_strand": run.segment.strand,
                        "a_maf_strand": a_row.strand,
                        "a_segment_id": "{}:{}-{}:{}".format(
                            a_row.seqid, run.segment.start0,
                            run.segment.end0, run.segment.strand),
                        "b_seqid": b_row.seqid,
                        "b_coord_min": _optional(b_evidence.coord_min),
                        "b_coord_max": _optional(b_evidence.coord_max),
                        "b_coord_first0": _optional(b_evidence.coord_first),
                        "b_coord_last0": _optional(b_evidence.coord_last),
                        "b_non_gap_bases": b_evidence.non_gap,
                        "b_cds_covered_bases": b_evidence.covered,
                        "b_strand_conflict_bases": b_evidence.conflict,
                        "b_cov_frac": "{:.6f}".format(b_cov),
                        "b_cds_strand": _evidence_strand(b_evidence),
                        "b_maf_strand": b_row.strand,
                        "c_seqid": c_row.seqid,
                        "c_coord_min": _optional(c_evidence.coord_min),
                        "c_coord_max": _optional(c_evidence.coord_max),
                        "c_coord_first0": _optional(c_evidence.coord_first),
                        "c_coord_last0": _optional(c_evidence.coord_last),
                        "c_non_gap_bases": c_evidence.non_gap,
                        "c_cds_covered_bases": c_evidence.covered,
                        "c_strand_conflict_bases": c_evidence.conflict,
                        "c_cov_frac": "{:.6f}".format(c_cov),
                        "c_cds_strand": _evidence_strand(c_evidence),
                        "c_maf_strand": c_row.strand,
                        "match_a": _format_match(primary.match_a),
                        "match_b": _format_match(primary.match_b),
                        "match_c": _format_match(primary.match_c),
                        "classification": primary.classification,
                        "classification_reason": primary.reason,
                    })
                    raw_rows_written += 1

    sensitivity = tuple(
        _summary_row(trio, threshold, counts[threshold],
                     blocks_seen, blocks_lt3)
        for threshold in thresholds
    )
    summary = next(row for row in sensitivity
                   if row["min_non_gap"] == min_non_gap)
    if summary["total_runs"] != raw_rows_written:
        raise RuntimeError(
            "{} summary/raw invariant failed: {} != {}".format(
                trio, summary["total_runs"], raw_rows_written)
        )
    return TrioAnalysis(summary=summary, sensitivity=sensitivity)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_input_manifest(data_root: Path, manifest_path: Path,
                            trios: Sequence[str]) -> List[Dict[str, Any]]:
    """Verify exact-assembly ingroup GFF files against the versioned manifest."""
    data_root = Path(data_root)
    selected = set(trios)
    validated = []
    roles_seen = {trio: set() for trio in selected}
    with Path(manifest_path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing_columns = set(MANIFEST_COLUMNS) - set(reader.fieldnames or ())
        if missing_columns:
            raise ValueError(
                "input manifest is missing columns: {}".format(
                    ", ".join(sorted(missing_columns)))
            )
        for row_number, row in enumerate(reader, start=2):
            trio = row["trio"]
            if trio not in selected:
                continue
            role = row["role"]
            if role not in ("ingroup B", "ingroup C"):
                raise ValueError(
                    "manifest row {} has unsupported role {!r}".format(
                        row_number, role))
            if role in roles_seen[trio]:
                raise ValueError(
                    "manifest contains duplicate {} {} entry".format(trio, role))
            roles_seen[trio].add(role)

            path = (data_root / "data" / "alignments" /
                    row["package_filename"])
            if not path.is_file():
                raise FileNotFoundError("manifest input is missing: {}".format(path))
            try:
                expected_bytes = int(row["bytes"])
            except ValueError as exc:
                raise ValueError(
                    "manifest row {} has invalid bytes".format(row_number)) from exc
            if path.stat().st_size != expected_bytes:
                raise ValueError(
                    "{} byte count mismatch: expected {}, observed {}".format(
                        path, expected_bytes, path.stat().st_size))
            observed_sha = _sha256(path)
            if observed_sha.lower() != row["sha256"].lower():
                raise ValueError(
                    "{} SHA-256 mismatch: expected {}, observed {}".format(
                        path, row["sha256"], observed_sha))

            accession_found = False
            cds_rows = 0
            seqids = set()
            try:
                opener = gzip.open if path.name.endswith(".gz") else open
                with opener(path, "rt") as gff_handle:
                    for line in gff_handle:
                        if line.startswith("#!genome-build-accession"):
                            accession_found = row["accession"] in line
                        if line.startswith("#"):
                            continue
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) == 9:
                            seqids.add(fields[0])
                            if fields[2] == "CDS":
                                cds_rows += 1
            except (OSError, EOFError) as exc:
                raise ValueError("invalid gzip/GFF input: {}".format(path)) from exc
            if not accession_found:
                raise ValueError(
                    "{} does not declare expected accession {}".format(
                        path, row["accession"]))
            if cds_rows == 0:
                raise ValueError("{} contains no CDS rows".format(path))
            checked = dict(row)
            checked.update({
                "path": path,
                "cds_rows": cds_rows,
                "seqids": frozenset(seqids),
            })
            validated.append(checked)

    for trio, roles in roles_seen.items():
        required = {"ingroup B", "ingroup C"}
        if roles != required:
            raise ValueError(
                "manifest must contain ingroup B and C for {}; observed {}"
                .format(trio, sorted(roles)))
    return validated


def _trio_paths(data_root: Path, trio: str) -> Tuple[Path, Path, Path, Path]:
    align = Path(data_root) / "data" / "alignments"
    return (
        align / "{}.maf.gz".format(trio),
        align / "{}_outgroup_CDS.gff.gz".format(trio),
        align / "{}_ingroupB_CDS.gff.gz".format(trio),
        align / "{}_ingroupC_CDS.gff.gz".format(trio),
    )


def _validate_seqid_compatibility(data_root: Path, trios: Sequence[str],
                                  records: Sequence[Dict[str, Any]]) -> None:
    by_key = {(row["trio"], row["role"]): row for row in records}
    for trio in trios:
        maf_path = _trio_paths(data_root, trio)[0]
        seqids = {"ingroup B": set(), "ingroup C": set()}
        with open_text(maf_path) as handle:
            for block_index, rows in enumerate(iter_maf_blocks(handle)):
                if len(rows) < 3:
                    continue
                if len(rows) != 3:
                    raise ValueError(
                        "MAF block {} must contain exactly 3 sequence rows; found {}"
                        .format(block_index, len(rows)))
                seqids["ingroup B"].add(rows[1].seqid)
                seqids["ingroup C"].add(rows[2].seqid)
        for role in ("ingroup B", "ingroup C"):
            missing = seqids[role] - set(by_key[(trio, role)]["seqids"])
            if missing:
                raise ValueError(
                    "{} {} MAF seqids are absent from GFF: {}".format(
                        trio, role, ", ".join(sorted(missing)[:10])))


def _write_table(path: Path, rows: Sequence[Dict[str, Any]]) -> None:
    with atomic_csv(path, SUMMARY_COLUMNS) as writer:
        for row in rows:
            writer.writerow(row)


def run_analysis(
    data_root: Path,
    trios: Sequence[str],
    out_dir: Path,
    min_non_gap: int = 10,
    sensitivity_thresholds: Sequence[int] = (1, 10, 30),
    manifest_path: Optional[Path] = None,
) -> Tuple[TrioAnalysis, ...]:
    """Validate inputs, run selected trios, and write combined result tables."""
    data_root = Path(data_root)
    out_dir = Path(out_dir)
    trios = tuple(trios)
    if not trios:
        raise ValueError("at least one trio is required")

    for trio in trios:
        for path in _trio_paths(data_root, trio):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(
                    "required input is missing or empty: {}".format(path))

    if manifest_path is not None:
        records = validate_input_manifest(data_root, manifest_path, trios)
        _validate_seqid_compatibility(data_root, trios, records)

    results = []
    for trio in trios:
        maf, a_gff, b_gff, c_gff = _trio_paths(data_root, trio)
        raw_path = out_dir / "cds_orth_strand_raw_{}.csv".format(trio)
        results.append(analyze_trio(
            trio,
            maf,
            a_gff,
            b_gff,
            c_gff,
            raw_path,
            min_non_gap=min_non_gap,
            sensitivity_thresholds=sensitivity_thresholds,
        ))

    summaries = [result.summary for result in results]
    sensitivity = [row for result in results for row in result.sensitivity]
    _write_table(out_dir / "cds_orth_strand_summary.csv", summaries)
    _write_table(out_dir / "cds_orth_strand_sensitivity.csv", sensitivity)
    return tuple(results)


def _parse_thresholds(value: str) -> Tuple[int, ...]:
    try:
        thresholds = tuple(int(item) for item in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "thresholds must be comma-separated integers") from exc
    if not thresholds or any(value < 1 for value in thresholds):
        raise argparse.ArgumentTypeError("thresholds must be positive")
    return thresholds


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", required=True, type=Path,
                        help="co-author package root containing data/alignments")
    parser.add_argument("--trios", nargs="+", default=list(DEFAULT_TRIOS),
                        help="trios to analyze")
    parser.add_argument("--min-nongap", type=int, default=10,
                        help="primary minimum B/C non-gap bases per run")
    parser.add_argument("--sensitivity-thresholds", type=_parse_thresholds,
                        default=(1, 10, 30),
                        help="comma-separated support thresholds")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="result table directory")
    parser.add_argument(
        "--input-manifest",
        type=Path,
        default=Path(__file__).resolve().parents[2] /
                "config" / "cds_orth_strand_inputs.tsv",
        help="versioned ingroup GFF provenance manifest",
    )
    args = parser.parse_args(argv)
    out_dir = (args.out_dir if args.out_dir is not None else
               args.data_root / "data" / "results" / "tables")
    results = run_analysis(
        args.data_root,
        args.trios,
        out_dir,
        min_non_gap=args.min_nongap,
        sensitivity_thresholds=args.sensitivity_thresholds,
        manifest_path=args.input_manifest,
    )
    writer = csv.DictWriter(sys.stdout, fieldnames=SUMMARY_COLUMNS)
    writer.writeheader()
    for result in results:
        writer.writerow(result.summary)


def iter_a_runs(block_index: int, coords: Sequence[Optional[int]],
                segments: Sequence[Any]) -> Iterator[AlignmentRun]:
    """Yield maximal consecutive non-gap columns in one A CDS segment."""
    if len(coords) != len(segments):
        raise ValueError("coordinate and segment arrays differ in length")
    active_segment = None
    start = None
    first_coord = None
    last_coord = None
    run_index = 0

    def make_run(end):
        return AlignmentRun(
            block_index=block_index,
            run_index=run_index,
            col_start0=start,
            col_end0=end,
            segment=active_segment,
            a_coord_first=first_coord,
            a_coord_last=last_coord,
        )

    for column, (coord, segment) in enumerate(zip(coords, segments)):
        if coord is None or segment is None:
            if active_segment is not None:
                yield make_run(column)
                run_index += 1
            active_segment = None
            start = first_coord = last_coord = None
            continue
        if active_segment is None:
            active_segment = segment
            start = column
            first_coord = last_coord = coord
            continue
        if segment != active_segment:
            yield make_run(column)
            run_index += 1
            active_segment = segment
            start = column
            first_coord = last_coord = coord
            continue
        last_coord = coord

    if active_segment is not None:
        yield make_run(len(coords))


@contextmanager
def atomic_csv(path: Path, columns: list) -> Iterator[csv.DictWriter]:
    """Write a CSV atomically, removing only the temporary file on failure."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        newline="",
        prefix=".{}.".format(path.name),
        suffix=".tmp",
        dir=str(path.parent),
        delete=False,
    )
    tmp_path = Path(handle.name)
    try:
        writer = csv.DictWriter(
            handle, fieldnames=columns, lineterminator="\n"
        )
        writer.writeheader()
        yield writer
        handle.flush()
        os.fsync(handle.fileno())
        handle.close()
        os.replace(str(tmp_path), str(path))
    except BaseException:
        handle.close()
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass
        raise


if __name__ == "__main__":
    main()
