import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path
from test_collect_run_summary import create_run
import collect_run_summary

ROOT = Path(__file__).resolve().parents[1]

class SpectrumShapeCliTest(unittest.TestCase):
    def test_multi_trio_outputs_and_exclusion(self):
        with tempfile.TemporaryDirectory(prefix="shape-", dir=ROOT / "tmp" / "codex") as work:
            work = Path(work)
            channels = [f"{a}[{s}]{b}" for s in ("C>A","C>G","C>T","T>A","T>C","T>G")
                        for a in "ACGT" for b in "ACGT"]
            for group in range(3):
                run = create_run(work / "input" / str(group))
                path = run / "metadata/metadata_manifest.json"
                manifest = json.loads(path.read_text())
                manifest["trio_key"] = f"trio-{group}"
                manifest["run_key"] = f"run-{group}"
                for artifact in manifest["artifacts"]:
                    if artifact["kind"] == "single":
                        content = "mutType\tmutNum\ttotalRootNum\n"
                        for i, channel in enumerate(channels):
                            count = (i + 1) if group == 0 else (97-i if group == 1 else (i % 16 + 1))
                            if group == 2 and artifact["slot"] == "org3" and i == 0:
                                count = 0
                            content += f"{channel}\t{count}\t10000\n"
                        table = run / artifact["path"]
                        table.write_text(content)
                        artifact["sha256"] = collect_run_summary.runs.sha256(table)
                path.write_text(json.dumps(manifest))
            command = ["Rscript", str(ROOT / "src/visualize/spectrum_shape.R"), str(work / "input"),
                       "--tsv-pattern", "*.tsv", "--taxonomy-levels", "phylum,class,order,species", "--output-dir", str(work / "output")]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            import csv
            def read(name):
                with (work / "output" / name).open() as handle:
                    return list(csv.DictReader(handle))
            self.assertEqual(len(read("clusters.csv")), 5)
            self.assertIn("taxonomy_species", read("samples.csv")[0])
            self.assertEqual(read("samples.csv")[0]["species"], "Homo beta")
            self.assertEqual(len(read("exclusions.csv")), 1)
            self.assertEqual(len(read("umap.csv")), 5)
            self.assertEqual(len(read("log2_matrix.csv")[0]), 97)
            for name in ("heatmap.pdf", "heatmap.png", "umap_phylum.png", "umap_cluster.pdf", "spectrum_cluster_1.png"):
                self.assertGreater((work / "output" / name).stat().st_size, 1000)
            first = (work / "output/umap.csv").read_text()
            for path in (work / "input").glob("*/*/metadata/metadata_manifest.json"):
                manifest = json.loads(path.read_text())
                for entry in manifest["organisms"]:
                    tax_path = path.parents[1] / entry["taxonomy_json"]
                    tax = json.loads(tax_path.read_text())
                    tax["reports"][0]["taxonomy"]["classification"]["phylum"] = {"name": "ChangedTaxon"}
                    tax_path.write_text(json.dumps(tax))
                    entry["taxonomy_json_sha256"] = collect_run_summary.runs.sha256(tax_path)
                path.write_text(json.dumps(manifest))
            command[-1] = str(work / "repeat")
            repeat = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(repeat.returncode, 0, repeat.stdout + repeat.stderr)
            self.assertEqual(first, (work / "repeat/umap.csv").read_text())
            if os.environ.get("SPECTRUM_SHAPE_TEST_ARTIFACTS"):
                import shutil
                shutil.copytree(work / "output", os.environ["SPECTRUM_SHAPE_TEST_ARTIFACTS"])
            for name in ("clusters.csv", "distances.csv", "log2_matrix.csv"):
                self.assertEqual((work / "output" / name).read_text(), (work / "repeat" / name).read_text())

if __name__ == "__main__":
    unittest.main()
