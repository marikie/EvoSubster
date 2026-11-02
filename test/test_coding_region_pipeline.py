import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("yq"), "pipeline requires yq")
class CodingRegionPipelineTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="coding pipeline ")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        (self.bin / "lastal").write_text("#!/bin/bash\necho stub LAST\n")
        (self.bin / "lastal").chmod(0o755)
        self.align = self.root / "align"
        self.align.mkdir()
        (self.align / "last_train.sh").write_text("#!/bin/bash\nexit 0\n")
        self.env = dict(os.environ, PATH=f"{self.bin}:{os.environ['PATH']}",
                        ALIGN_DIR_OVERRIDE=str(self.align))
        self.fasta = []
        self.gffs = []
        for name in ("Alpha_one", "Beta_two", "Gamma_three"):
            genome = self.root / name
            genome.mkdir()
            fasta = genome / "genomic.fna"
            fasta.write_text(">chr1\nACGTNNNN\n")
            self.fasta.append(str(fasta))
            self.gffs.append(genome / "genomic.gff")
        self.out = self.root / "out"
        self.misc = self.out / "Alpone1_Bettwo2_Gamthr3/20260101/statistics/misc"

    def annotate(self, i):
        self.gffs[i].write_text("chr1\ttest\tCDS\t1\t2\t.\t+\t0\tID=cds\n")

    def run_pipeline(self, *options, outgroup=None):
        return subprocess.run(
            ["bash", str(ROOT / "src/sbst.sh"), "20260101", *self.fasta,
             str(outgroup) if outgroup is not None else "NO_GFF_FILE",
             "--out-dir", str(self.out), "--idt-only", *options],
            env=self.env, capture_output=True, text=True,
        )

    def test_all_species_and_reuse_of_completed_results(self):
        for i in range(3):
            self.annotate(i)
        result = self.run_pipeline(outgroup=self.gffs[0])
        self.assertEqual(result.returncode, 0, result.stderr)
        files = sorted(self.misc.glob("*_codingRegionPercentage_*.out"))
        self.assertEqual(len(files), 3, result.stdout)
        for path in files:
            self.assertEqual(path.read_text(), "Coding region percentage: 25.00%\n")
        stamps = [p.stat().st_mtime_ns for p in files]
        result = self.run_pipeline(outgroup=self.gffs[0])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(stamps, [p.stat().st_mtime_ns for p in files])

    def test_missing_annotations_and_outgroup_sentinel(self):
        self.annotate(0)
        self.annotate(2)
        result = self.run_pipeline()
        self.assertEqual(result.returncode, 0, result.stderr)
        files = list(self.misc.glob("*_codingRegionPercentage_*.out"))
        self.assertEqual([p.name for p in files],
                         ["Gamthr3_codingRegionPercentage_20260101.out"])

    def test_explicit_ingroup_gff_paths_and_empty_result_recovery(self):
        annotation = self.root / "custom annotation.gff"
        annotation.write_text("chr1\ttest\tCDS\t1\t4\t.\t+\t0\tID=c\n")
        self.misc.mkdir(parents=True)
        output = self.misc / "Bettwo2_codingRegionPercentage_20260101.out"
        output.touch()
        result = self.run_pipeline("--org2-gff", str(annotation),
                                   f"--org3-gff={annotation}")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(output.read_text(), "Coding region percentage: 50.00%\n")
        self.assertEqual(len(list(self.misc.glob("*_codingRegionPercentage_*.out"))), 2)

    def test_invalid_annotation_stops_pipeline_without_result(self):
        self.gffs[1].write_text("missing\ttest\tCDS\t1\t2\t.\t+\t0\tID=c\n")
        result = self.run_pipeline()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(list(self.misc.glob("*_codingRegionPercentage_*.out")))

    def test_missing_explicit_annotation_is_error(self):
        result = self.run_pipeline("--org2-gff", str(self.root / "absent.gff"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("GFF", result.stderr)

    def test_empty_outgroup_argument_is_not_autodetected(self):
        self.annotate(0)
        result = self.run_pipeline(outgroup="")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("GFF", result.stderr)

    def test_empty_option_is_error(self):
        for option in ("--org2-gff=", "--org3-gff="):
            with self.subTest(option=option):
                result = self.run_pipeline(option)
                self.assertNotEqual(result.returncode, 0)

    def test_explicit_sentinel_disables_ingroup_autodetection(self):
        self.annotate(1)
        result = self.run_pipeline("--org2-gff", "NO_GFF_FILE")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(list(self.misc.glob("*_codingRegionPercentage_*.out")))


if __name__ == "__main__":
    unittest.main()
