import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src/metrics/coding_region_percentage.py"


def cds(seq, start, end, strand="+"):
    return f"{seq}\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds\n"


class CodingRegionPercentageTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="coding region ")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fasta = self.root / "genome.fa"
        self.gff = self.root / "genomic.gff"
        self.output = self.root / "percentage.out"

    def run_metric(self, fasta, gff):
        self.fasta.write_text(fasta)
        self.gff.write_text(gff)
        return subprocess.run(
            [sys.executable, str(SCRIPT), str(self.fasta), str(self.gff),
             "--output", str(self.output)],
            capture_output=True, text=True,
        )

    def test_union_counts_overlaps_duplicates_strands_and_unannotated_contigs(self):
        result = self.run_metric(
            ">chr1 description\nACGTN\nNNNNN\n>chr2\nACGTN\n>unannotated\nNNNNN\n",
            "##gff-version 3\n" + cds("chr1", 2, 5) + cds("chr1", 4, 7, "-")
            + cds("chr1", 2, 5) + cds("chr2", 5, 5),
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.output.read_text(), "Coding region percentage: 35.00%\n")

    def test_no_cds_is_zero_and_embedded_fasta_is_ignored(self):
        result = self.run_metric(
            ">chr1\nACGT\n",
            "##gff-version 3\nchr1\ttest\tgene\t1\t4\t.\t+\t.\tID=gene\n"
            "##FASTA\n>chr1\nACGT\n",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.output.read_text(), "Coding region percentage: 0.00%\n")

    def test_adjacent_cds_cover_entire_genome_with_inclusive_boundaries(self):
        result = self.run_metric(">chr1\nACGT\n", cds("chr1", 1, 2) + cds("chr1", 3, 4))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(self.output.read_text(), "Coding region percentage: 100.00%\n")

    def test_invalid_cds_fails_without_output(self):
        for annotation in (cds("absent", 1, 2), cds("chr1", 0, 2),
                           cds("chr1", 3, 2), cds("chr1", 1, 5),
                           cds("chr1", "bad", 2), "chr1\ttest\tCDS\t1\n"):
            with self.subTest(annotation=annotation):
                result = self.run_metric(">chr1\nACGT\n", annotation)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Error:", result.stderr)
                self.assertFalse(self.output.exists())

    def test_invalid_fasta_fails(self):
        for fasta in ("", ">chr1\n", "ACGT\n", ">chr1\nAC\n>chr1\nGT\n"):
            with self.subTest(fasta=fasta):
                result = self.run_metric(fasta, cds("chr1", 1, 2))
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("Error:", result.stderr)
                self.assertFalse(self.output.exists())

    def test_failed_recalculation_preserves_existing_output(self):
        self.output.write_text("previous result\n")
        result = self.run_metric(">chr1\nACGT\n", cds("chr1", 1, 5))
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.output.read_text(), "previous result\n")
        self.assertEqual(sorted(p.name for p in self.root.iterdir()),
                         ["genome.fa", "genomic.gff", "percentage.out"])

    def test_output_replace_failure_cleans_staging_file(self):
        from unittest.mock import patch

        spec = importlib.util.spec_from_file_location("coding_region_percentage", SCRIPT)
        metric = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(metric)
        self.fasta.write_text(">chr1\nACGT\n")
        self.gff.write_text(cds("chr1", 1, 2))
        with patch.object(sys, "argv", [str(SCRIPT), str(self.fasta), str(self.gff),
                                       "--output", str(self.output)]), \
                patch.object(metric.os, "replace", side_effect=OSError("replace failed")):
            self.assertEqual(metric.main(), 1)
        self.assertFalse(self.output.exists())
        self.assertEqual(len(list(self.root.iterdir())), 2)


if __name__ == "__main__":
    unittest.main()
