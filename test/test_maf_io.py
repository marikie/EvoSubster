#!/usr/bin/env python3
"""Tests for strict shared MAF parsing and coordinate helpers."""

import gzip
import os
import sys
import tempfile
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from src.align.maf_io import (
    MafRow,
    aligned_forward_coords,
    forward_coord,
    iter_maf_blocks,
    open_text,
)


def maf_row(
    seqid="ref",
    start0=10,
    span=4,
    strand="+",
    seqlen=100,
    sequence="ACGT",
):
    return "s {} {} {} {} {} {}\n".format(
        seqid, start0, span, strand, seqlen, sequence
    )


class TestOpenText(unittest.TestCase):
    def test_reads_plain_and_gzip_text(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            plain_path = os.path.join(temp_dir, "input.maf")
            gzip_path = plain_path + ".gz"
            with open(plain_path, "w") as handle:
                handle.write("plain\n")
            with gzip.open(gzip_path, "wt") as handle:
                handle.write("compressed\n")

            with open_text(plain_path) as handle:
                self.assertEqual(handle.read(), "plain\n")
            with open_text(gzip_path) as handle:
                self.assertEqual(handle.read(), "compressed\n")


class TestIterMafBlocks(unittest.TestCase):
    def test_preserves_all_row_fields_and_block_boundaries(self):
        lines = [
            "##maf version=1\n",
            "a score=1\n",
            maf_row("ref", 10, 3, "+", 100, "AC-G"),
            maf_row("query", 20, 3, "-", 200, "A-TG"),
            "\n",
            "a score=2\n",
            maf_row("other", 5, 2, "+", 50, "TT"),
        ]

        blocks = list(iter_maf_blocks(lines))

        self.assertEqual(
            blocks,
            [
                [
                    MafRow("ref", 10, 3, "+", 100, "AC-G"),
                    MafRow("query", 20, 3, "-", 200, "A-TG"),
                ],
                [MafRow("other", 5, 2, "+", 50, "TT")],
            ],
        )

    def test_new_alignment_header_ends_previous_block_without_blank_line(self):
        lines = [
            "a score=1\n",
            maf_row("first", sequence="ACGT"),
            "a score=2\n",
            maf_row("second", sequence="TGCA"),
        ]

        blocks = list(iter_maf_blocks(lines))

        self.assertEqual([[row.seqid for row in block] for block in blocks], [
            ["first"],
            ["second"],
        ])

    def test_rejects_malformed_sequence_row(self):
        with self.assertRaisesRegex(ValueError, "expected 7 fields"):
            list(iter_maf_blocks(["a\n", "s ref 0 4 + 100\n"]))

    def test_rejects_invalid_strand(self):
        with self.assertRaisesRegex(ValueError, "strand"):
            list(iter_maf_blocks(["a\n", maf_row(strand="?")]))

    def test_rejects_coordinates_outside_source_sequence(self):
        with self.assertRaisesRegex(ValueError, "source sequence"):
            list(iter_maf_blocks(["a\n", maf_row(start0=98, span=4)]))

    def test_rejects_span_mismatch(self):
        with self.assertRaisesRegex(ValueError, "span"):
            list(iter_maf_blocks(["a\n", maf_row(span=4, sequence="AC-G")]))

    def test_rejects_unequal_aligned_lengths(self):
        lines = [
            "a\n",
            maf_row("ref", span=4, sequence="ACGT"),
            maf_row("query", span=3, sequence="ACG"),
            "\n",
        ]
        with self.assertRaisesRegex(ValueError, "aligned length"):
            list(iter_maf_blocks(lines))


class TestCoordinates(unittest.TestCase):
    def test_forward_coord_on_plus_strand(self):
        row = MafRow("ref", 10, 4, "+", 100, "ACGT")
        self.assertEqual(forward_coord(row, 0), 10)
        self.assertEqual(forward_coord(row, 3), 13)

    def test_forward_coord_on_minus_strand(self):
        row = MafRow("query", 10, 4, "-", 100, "ACGT")
        self.assertEqual(forward_coord(row, 0), 89)
        self.assertEqual(forward_coord(row, 3), 86)

    def test_rejects_ungapped_offset_outside_row_span(self):
        row = MafRow("ref", 10, 4, "+", 100, "ACGT")
        with self.assertRaisesRegex(ValueError, "ungapped offset"):
            forward_coord(row, 4)

    def test_aligned_forward_coords_include_gaps(self):
        plus = MafRow("ref", 10, 3, "+", 100, "A-CG")
        minus = MafRow("query", 10, 3, "-", 100, "A-CG")

        self.assertEqual(aligned_forward_coords(plus), [10, None, 11, 12])
        self.assertEqual(aligned_forward_coords(minus), [89, None, 88, 87])


if __name__ == "__main__":
    unittest.main(verbosity=2)
