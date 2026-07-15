#!/usr/bin/env python3
"""Tests for the independent CDS/strand consistency QC."""

import csv
import gzip
import hashlib
import os
import subprocess
import sys
import tempfile
import unittest
from dataclasses import dataclass
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPO = HERE.parent
ANALYSIS = REPO / "src" / "analysis"
sys.path.insert(0, str(ANALYSIS))

import cds_orth_strand_check as C  # noqa: E402


@dataclass(frozen=True)
class FakeSegment:
    start0: int
    end0: int
    strand: str


class FakeIndex:
    def __init__(self, segments):
        self._segments = tuple(segments)

    def segments_for(self, seqid):
        return self._segments if seqid == "seq" else ()


def evidence(**overrides):
    values = {
        "a_cds_strand": "+",
        "a_maf_strand": "+",
        "b_maf_strand": "+",
        "c_maf_strand": "+",
        "b_non_gap": 12,
        "b_covered": 12,
        "b_plus": 12,
        "b_minus": 0,
        "b_conflict": 0,
        "c_non_gap": 12,
        "c_covered": 12,
        "c_plus": 12,
        "c_minus": 0,
        "c_conflict": 0,
    }
    values.update(overrides)
    return values


class TestClassification(unittest.TestCase):
    def test_consistent_when_all_matches_true(self):
        result = C.classify_evidence(**evidence(), min_non_gap=10)
        self.assertEqual(result.classification, C.CLASS_CONSISTENT)
        self.assertEqual((result.match_a, result.match_b, result.match_c),
                         (True, True, True))

    def test_consistent_when_all_matches_false(self):
        result = C.classify_evidence(**evidence(
            a_cds_strand="-",
            b_minus=12,
            b_plus=0,
            c_minus=12,
            c_plus=0,
        ), min_non_gap=10)
        self.assertEqual(result.classification, C.CLASS_CONSISTENT)
        self.assertEqual((result.match_a, result.match_b, result.match_c),
                         (False, False, False))

    def test_low_cds_coverage_is_noncoding(self):
        result = C.classify_evidence(**evidence(b_covered=5, b_plus=5),
                                     min_non_gap=10)
        self.assertEqual(result.classification, C.CLASS_NONCODING)
        self.assertEqual(result.reason, "b_or_c_below_cds_coverage")

    def test_mixed_strand_matches_are_contradictory(self):
        result = C.classify_evidence(**evidence(b_plus=0, b_minus=12),
                                     min_non_gap=10)
        self.assertEqual(result.classification, C.CLASS_CONTRADICTION)

    def test_low_non_gap_support_is_insufficient(self):
        result = C.classify_evidence(**evidence(b_non_gap=9, b_covered=9,
                                                b_plus=9), min_non_gap=10)
        self.assertEqual(result.classification, C.CLASS_INSUFFICIENT)
        self.assertEqual(result.reason, "b_or_c_below_min_non_gap")

    def test_opposite_strand_overlap_is_insufficient(self):
        result = C.classify_evidence(**evidence(c_plus=11, c_conflict=1),
                                     min_non_gap=10)
        self.assertEqual(result.classification, C.CLASS_INSUFFICIENT)
        self.assertEqual(result.reason, "b_or_c_cds_strand_conflict")

    def test_tied_strand_vote_is_insufficient(self):
        result = C.classify_evidence(**evidence(c_plus=6, c_minus=6),
                                     min_non_gap=10)
        self.assertEqual(result.classification, C.CLASS_INSUFFICIENT)
        self.assertEqual(result.reason, "b_or_c_cds_strand_tie")

    def test_ambiguous_a_segment_is_insufficient(self):
        result = C.classify_evidence(**evidence(a_cds_strand="."),
                                     min_non_gap=10)
        self.assertEqual(result.classification, C.CLASS_INSUFFICIENT)
        self.assertEqual(result.reason, "a_cds_strand_ambiguous")


class TestAtomicCsv(unittest.TestCase):
    def test_atomic_csv_removes_partial_target_on_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "rows.csv"
            with self.assertRaisesRegex(RuntimeError, "stop"):
                with C.atomic_csv(target, ["value"]) as writer:
                    writer.writerow({"value": 1})
                    raise RuntimeError("stop")
            self.assertFalse(target.exists())
            self.assertEqual(list(Path(tmp).glob(".rows.csv.*.tmp")), [])

    def test_atomic_csv_replaces_target_after_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "rows.csv"
            target.write_text("old\n")
            with C.atomic_csv(target, ["value"]) as writer:
                writer.writerow({"value": 7})
            with target.open(newline="") as fh:
                self.assertEqual(list(csv.DictReader(fh)), [{"value": "7"}])

    def test_atomic_csv_uses_lf_line_endings(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "rows.csv"
            with C.atomic_csv(target, ["value"]) as writer:
                writer.writerow({"value": 7})
            self.assertEqual(target.read_bytes(), b"value\n7\n")


class TestMonotonicAnnotation(unittest.TestCase):
    def setUp(self):
        self.plus = FakeSegment(0, 3, "+")
        self.conflict = FakeSegment(3, 5, ".")
        self.minus = FakeSegment(5, 8, "-")
        self.index = FakeIndex([self.plus, self.conflict, self.minus])

    def test_annotates_increasing_coordinates_with_gaps(self):
        states = C.annotate_monotonic(
            self.index, "seq", [1, 2, None, 3, 4, 5, 6]
        )
        self.assertEqual(states, [
            self.plus, self.plus, None, self.conflict, self.conflict,
            self.minus, self.minus,
        ])

    def test_annotates_decreasing_coordinates(self):
        states = C.annotate_monotonic(
            self.index, "seq", [6, 5, 4, 3, 2, 1]
        )
        self.assertEqual(states, [
            self.minus, self.minus, self.conflict, self.conflict,
            self.plus, self.plus,
        ])

    def test_collects_species_evidence_and_coordinate_bounds(self):
        result = C.collect_species_evidence(
            self.index, "seq", [1, 2, None, 3, 4, 5, 6]
        )
        self.assertEqual(
            (result.non_gap, result.covered, result.plus, result.minus,
             result.conflict),
            (6, 6, 2, 2, 2),
        )
        self.assertEqual(
            (result.coord_min, result.coord_max, result.coord_first,
             result.coord_last),
            (1, 6, 1, 6),
        )

    def test_a_gap_breaks_a_run_even_with_same_segment(self):
        coords = [0, 1, None, 2, 3, 5]
        states = C.annotate_monotonic(self.index, "seq", coords)
        runs = list(C.iter_a_runs(7, coords, states))
        self.assertEqual(
            [(r.run_index, r.col_start0, r.col_end0, r.segment.strand)
             for r in runs],
            [(0, 0, 2, "+"), (1, 3, 4, "+"), (2, 4, 5, "."),
             (3, 5, 6, "-")],
        )


def write_text(path, text):
    path.write_text(text)
    return path


def write_gzip_text(path, text):
    with gzip.open(path, "wt") as fh:
        fh.write(text)
    return path


def maf_block(a_start, a_seq, b_start, b_seq, c_start, c_seq,
              b_strand="+", c_strand="+"):
    return "\n".join([
        "a",
        "s A {} {} + 1000 {}".format(
            a_start, len(a_seq.replace("-", "")), a_seq),
        "s B {} {} {} 1000 {}".format(
            b_start, len(b_seq.replace("-", "")), b_strand, b_seq),
        "s C {} {} {} 1000 {}".format(
            c_start, len(c_seq.replace("-", "")), c_strand, c_seq),
        "",
    ])


def build_single_run_package(root):
    align = root / "data" / "alignments"
    align.mkdir(parents=True)
    twelve = "A" * 12
    write_gzip_text(
        align / "toy.maf.gz",
        maf_block(0, twelve, 0, twelve, 0, twelve),
    )
    write_gzip_text(
        align / "toy_outgroup_CDS.gff.gz",
        "A\ttest\tCDS\t1\t12\t.\t+\t0\tParent=a\n",
    )
    records = []
    for role, seqid, accession in (
        ("B", "B", "GCA_000000001.1"),
        ("C", "C", "GCA_000000002.1"),
    ):
        path = align / "toy_ingroup{}_CDS.gff.gz".format(role)
        write_gzip_text(
            path,
            "#!genome-build-accession NCBI_Assembly:{}\n".format(accession)
            + "{}\ttest\tCDS\t1\t12\t.\t+\t0\tParent=x\n".format(seqid),
        )
        records.append((role, accession, path,
                        hashlib.sha256(path.read_bytes()).hexdigest()))
    manifest = root / "inputs.tsv"
    lines = ["\t".join([
        "trio", "role", "accession", "annotation_source",
        "package_filename", "bytes", "sha256", "url",
    ])]
    for role, accession, path, digest in records:
        lines.append("\t".join([
            "toy", "ingroup " + role, accession, "test", path.name,
            str(path.stat().st_size), digest,
            "https://example.invalid/toy.gff.gz",
        ]))
    manifest.write_text("\n".join(lines) + "\n")
    return manifest


class TestSyntheticTrio(unittest.TestCase):
    def test_streams_four_classes_and_sensitivity_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            twelve = "A" * 12
            maf = write_text(root / "toy.maf", "\n".join([
                maf_block(0, twelve, 0, twelve, 0, twelve),
                maf_block(20, twelve, 20, twelve, 20, twelve),
                maf_block(40, twelve, 40, twelve, 40, twelve),
                maf_block(60, twelve, 60, "-" * 9 + "AAA", 60, twelve),
                "a\ns A 80 3 + 1000 AAA\ns B 80 3 + 1000 AAA\n",
            ]))
            a_gff = write_text(root / "a.gff", "\n".join([
                "A\ttest\tCDS\t1\t12\t.\t+\t0\tParent=a1",
                "A\ttest\tCDS\t21\t32\t.\t+\t0\tParent=a2",
                "A\ttest\tCDS\t41\t52\t.\t+\t0\tParent=a3",
                "A\ttest\tCDS\t61\t72\t.\t+\t0\tParent=a4",
                "",
            ]))
            b_gff = write_text(root / "b.gff", "\n".join([
                "B\ttest\tCDS\t1\t12\t.\t+\t0\tParent=b1",
                "B\ttest\tCDS\t41\t52\t.\t-\t0\tParent=b3",
                "B\ttest\tCDS\t61\t63\t.\t+\t0\tParent=b4",
                "",
            ]))
            c_gff = write_text(root / "c.gff", "\n".join([
                "C\ttest\tCDS\t1\t12\t.\t+\t0\tParent=c1",
                "C\ttest\tCDS\t21\t32\t.\t+\t0\tParent=c2",
                "C\ttest\tCDS\t41\t52\t.\t+\t0\tParent=c3",
                "C\ttest\tCDS\t61\t72\t.\t+\t0\tParent=c4",
                "",
            ]))
            raw_path = root / "raw.csv"

            result = C.analyze_trio(
                "toy", maf, a_gff, b_gff, c_gff, raw_path,
                min_non_gap=10,
                sensitivity_thresholds=(1, 10, 30),
            )

            with raw_path.open(newline="") as fh:
                raw = list(csv.DictReader(fh))

        self.assertEqual([row["classification"] for row in raw], [
            C.CLASS_CONSISTENT,
            C.CLASS_NONCODING,
            C.CLASS_CONTRADICTION,
            C.CLASS_INSUFFICIENT,
        ])
        self.assertEqual(len(raw), 4)
        self.assertEqual(result.summary["total_runs"], 4)
        self.assertEqual(result.summary["classified_runs"], 3)
        self.assertEqual(result.summary["blocks_seen"], 5)
        self.assertEqual(result.summary["blocks_lt3"], 1)
        by_threshold = {row["min_non_gap"]: row for row in result.sensitivity}
        self.assertEqual(by_threshold[1][C.CLASS_CONSISTENT], 2)
        self.assertEqual(by_threshold[10][C.CLASS_INSUFFICIENT], 1)
        self.assertEqual(by_threshold[30][C.CLASS_INSUFFICIENT], 4)

    def test_four_row_block_fails_without_partial_raw_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            maf = write_text(root / "bad.maf", "\n".join([
                "a",
                "s A 0 3 + 100 AAA",
                "s B 0 3 + 100 AAA",
                "s C 0 3 + 100 AAA",
                "s D 0 3 + 100 AAA",
                "",
            ]))
            gff = write_text(
                root / "x.gff",
                "A\ttest\tCDS\t1\t3\t.\t+\t0\tParent=x\n",
            )
            raw_path = root / "raw.csv"
            with self.assertRaisesRegex(ValueError, "exactly 3"):
                C.analyze_trio("bad", maf, gff, gff, gff, raw_path)
            self.assertFalse(raw_path.exists())


class TestManifestAndRun(unittest.TestCase):
    def test_validates_manifest_checksum_accession_and_cds(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            align = root / "data" / "alignments"
            align.mkdir(parents=True)
            records = []
            for role, seqid, accession in (
                ("B", "B", "GCA_000000001.1"),
                ("C", "C", "GCA_000000002.1"),
            ):
                gff = align / "toy_ingroup{}_CDS.gff.gz".format(role)
                with gzip.open(gff, "wt") as fh:
                    fh.write("#!genome-build-accession NCBI_Assembly:{}\n".format(
                        accession))
                    fh.write("{}\ttest\tCDS\t1\t3\t.\t+\t0\tParent=x\n".format(
                        seqid))
                digest = hashlib.sha256(gff.read_bytes()).hexdigest()
                records.append((role, accession, gff, digest))
            manifest = root / "inputs.tsv"
            header = "\t".join([
                "trio", "role", "accession", "annotation_source",
                "package_filename", "bytes", "sha256", "url",
            ])
            body = []
            for role, accession, gff, digest in records:
                body.append("\t".join([
                    "toy", "ingroup " + role, accession, "test",
                    gff.name, str(gff.stat().st_size), digest,
                    "https://example.invalid/toy.gff.gz",
                ]))
            manifest.write_text(header + "\n" + "\n".join(body) + "\n")

            validated = C.validate_input_manifest(root, manifest, ("toy",))

            self.assertEqual(len(validated), 2)
            self.assertEqual([record["cds_rows"] for record in validated], [1, 1])

            text = manifest.read_text().replace(records[0][3], "0" * 64)
            manifest.write_text(text)
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                C.validate_input_manifest(root, manifest, ("toy",))

    def test_run_analysis_writes_combined_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = build_single_run_package(root)
            out_dir = root / "results"

            results = C.run_analysis(
                root,
                ("toy",),
                out_dir,
                min_non_gap=10,
                sensitivity_thresholds=(1, 10, 30),
                manifest_path=manifest,
            )

            summary_path = out_dir / "cds_orth_strand_summary.csv"
            sensitivity_path = out_dir / "cds_orth_strand_sensitivity.csv"
            self.assertTrue(summary_path.is_file())
            self.assertTrue(sensitivity_path.is_file())
            with summary_path.open(newline="") as fh:
                summary = list(csv.DictReader(fh))
            with sensitivity_path.open(newline="") as fh:
                sensitivity = list(csv.DictReader(fh))
            self.assertEqual(len(results), 1)
            self.assertEqual(summary[0][C.CLASS_CONSISTENT], "1")
            self.assertEqual([row["min_non_gap"] for row in sensitivity],
                             ["1", "10", "30"])

    def test_cli_executes_after_all_helpers_are_defined(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = build_single_run_package(root)
            out_dir = root / "results"
            script = REPO / "src" / "analysis" / "cds_orth_strand_check.py"

            completed = subprocess.run([
                sys.executable,
                str(script),
                "--data-root", str(root),
                "--trios", "toy",
                "--out-dir", str(out_dir),
                "--input-manifest", str(manifest),
            ], cwd=str(REPO), capture_output=True, text=True)

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue((out_dir / "cds_orth_strand_summary.csv").is_file())


if __name__ == "__main__":
    unittest.main(verbosity=2)
