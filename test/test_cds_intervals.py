#!/usr/bin/env python3
"""Tests for deterministic CDS interval segmentation."""

import gzip
import os
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

try:
    from src.analysis.cds_intervals import CDSSegment, CDSIndex
except ModuleNotFoundError as exc:
    CDSSegment = None
    CDSIndex = None
    _IMPORT_ERROR = exc
else:
    _IMPORT_ERROR = None


def gff_record(seqid, feature_type, start1, end1, strand, attributes="ID=x"):
    """Return one GFF3 record."""
    fields = (
        seqid,
        "test",
        feature_type,
        str(start1),
        str(end1),
        ".",
        strand,
        "0",
        attributes,
    )
    return "\t".join(fields) + "\n"


class TestCDSModuleAvailable(unittest.TestCase):
    def test_cds_interval_module_is_available(self):
        self.assertIsNone(_IMPORT_ERROR)


@unittest.skipIf(_IMPORT_ERROR is not None, "production module not implemented")
class TestCDSIndex(unittest.TestCase):
    def _write_gff(self, contents, suffix=".gff3"):
        handle = tempfile.NamedTemporaryFile(
            mode="w", suffix=suffix, delete=False, encoding="utf-8"
        )
        with handle:
            handle.write(contents)
        self.addCleanup(os.unlink, handle.name)
        return handle.name

    def _write_gzip_gff(self, contents):
        handle = tempfile.NamedTemporaryFile(suffix=".gff3.gz", delete=False)
        handle.close()
        with gzip.open(handle.name, "wt", encoding="utf-8") as stream:
            stream.write(contents)
        self.addCleanup(os.unlink, handle.name)
        return handle.name

    def test_converts_gff_coordinates_and_exposes_sorted_seqids(self):
        path = self._write_gff(
            "##gff-version 3\n"
            + gff_record("chr2", "CDS", 3, 5, "+")
            + gff_record("chr1", "CDS", 10, 12, "-")
        )

        index = CDSIndex.from_gff(path)

        self.assertEqual(index.seqids, ("chr1", "chr2"))
        self.assertEqual(
            index.segments_for("chr2"), (CDSSegment(2, 5, "+"),)
        )
        self.assertIsNone(index.query("chr2", 1))
        self.assertEqual(index.query("chr2", 2), CDSSegment(2, 5, "+"))
        self.assertEqual(index.query("chr2", 4), CDSSegment(2, 5, "+"))
        self.assertIsNone(index.query("chr2", 5))

    def test_merges_touching_same_strand_intervals_but_preserves_gaps(self):
        path = self._write_gff(
            gff_record("chr1", "CDS", 1, 3, "+", "Parent=tx1")
            + gff_record("chr1", "CDS", 4, 6, "+", "Parent=tx2")
            + gff_record("chr1", "CDS", 8, 9, "+", "Parent=tx3")
        )

        index = CDSIndex.from_gff(path)

        self.assertEqual(
            index.segments_for("chr1"),
            (CDSSegment(0, 6, "+"), CDSSegment(7, 9, "+")),
        )
        self.assertIsNone(index.query("chr1", 6))

    def test_touching_opposite_strands_are_not_ambiguous(self):
        path = self._write_gff(
            gff_record("chr1", "CDS", 1, 3, "+")
            + gff_record("chr1", "CDS", 4, 6, "-")
        )

        index = CDSIndex.from_gff(path)

        self.assertEqual(
            index.segments_for("chr1"),
            (CDSSegment(0, 3, "+"), CDSSegment(3, 6, "-")),
        )

    def test_unions_same_strand_overlap_and_duplicate_parent_records(self):
        duplicate = gff_record("chr1", "CDS", 1, 5, "+", "Parent=tx1")
        path = self._write_gff(
            duplicate
            + duplicate
            + gff_record("chr1", "CDS", 4, 8, "+", "Parent=tx2")
        )

        index = CDSIndex.from_gff(path)

        self.assertEqual(
            index.segments_for("chr1"), (CDSSegment(0, 8, "+"),)
        )

    def test_marks_only_opposite_strand_overlap_as_ambiguous(self):
        path = self._write_gff(
            gff_record("chr1", "CDS", 1, 6, "+")
            + gff_record("chr1", "CDS", 4, 9, "-")
        )

        index = CDSIndex.from_gff(path)

        self.assertEqual(
            index.segments_for("chr1"),
            (
                CDSSegment(0, 3, "+"),
                CDSSegment(3, 6, "."),
                CDSSegment(6, 9, "-"),
            ),
        )
        self.assertEqual(index.query("chr1", 3).strand, ".")
        self.assertEqual(index.query("chr1", 5).strand, ".")
        self.assertEqual(index.query("chr1", 6).strand, "-")

    def test_absent_seqid_returns_empty_segments_and_no_query_match(self):
        path = self._write_gff(gff_record("chr1", "CDS", 1, 3, "+"))
        index = CDSIndex.from_gff(path)

        self.assertEqual(index.segments_for("missing"), ())
        self.assertIsNone(index.query("missing", 0))

    def test_reads_gzip_and_honors_custom_feature_type(self):
        path = self._write_gzip_gff(
            gff_record("chr1", "CDS", 1, 3, "+")
            + gff_record("chr1", "exon", 5, 7, "-")
        )

        index = CDSIndex.from_gff(path, feature_type="exon")

        self.assertEqual(
            index.segments_for("chr1"), (CDSSegment(4, 7, "-"),)
        )

    def test_rejects_invalid_selected_feature_coordinates(self):
        invalid_records = (
            gff_record("chr1", "CDS", 0, 3, "+"),
            gff_record("chr1", "CDS", 5, 4, "+"),
            gff_record("chr1", "CDS", "bad", 4, "+"),
        )
        for record in invalid_records:
            with self.subTest(record=record):
                path = self._write_gff(record)
                with self.assertRaisesRegex(ValueError, "line 1"):
                    CDSIndex.from_gff(path)

    def test_rejects_invalid_selected_feature_strand(self):
        for strand in (".", "?"):
            with self.subTest(strand=strand):
                path = self._write_gff(
                    gff_record("chr1", "CDS", 1, 3, strand)
                )
                with self.assertRaisesRegex(ValueError, "line 1"):
                    CDSIndex.from_gff(path)

    def test_rejects_negative_query_coordinates(self):
        path = self._write_gff(gff_record("chr1", "CDS", 1, 3, "+"))
        index = CDSIndex.from_gff(path)

        with self.assertRaises(ValueError):
            index.query("chr1", -1)


if __name__ == "__main__":
    unittest.main()
