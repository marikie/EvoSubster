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
    def test_mixed_v2_alias_and_timestamp(self):
        with tempfile.TemporaryDirectory(prefix="mixed-", dir=ROOT / "tmp/codex") as work:
            work = Path(work)
            run = create_run(work / "v2")
            manifest = json.loads((run / "metadata/metadata_manifest.json").read_text())
            trio = work / "cohort/trio"
            old = trio / "bench_20260101/metadata"
            old.mkdir(parents=True)
            (old / "metadata_manifest.json").write_text(json.dumps({"organisms": manifest["organisms"]}))
            legacy_counts = old.parent / "statistics/singlenuc"
            legacy_counts.mkdir(parents=True)
            channels_before = [f"{a}[{v}]{b}" for v in ("C>A","C>G","C>T","T>A","T>C","T>G") for a in "ACGT" for b in "ACGT"]
            for entry in manifest["organisms"]:
                if entry["role"] == "ingroup":
                    name = entry["accession"] + "_" + entry["short_name"] + "_old_ncds.tsv"
                    (legacy_counts / name).write_text("mutType\tmutNum\ttotalRootNum\n" + "".join(f"{c}\t1\t100\n" for c in channels_before))
            alias = trio / "bench_20260101_120000"
            alias.symlink_to(run.resolve(), target_is_directory=True)
            expression = "source('src/visualize/spectrum_shape.R'); x<-shape_mixed_ingest(commandArgs(TRUE)[1], 'whole', 'phylum', 'ward', 'bench_*', file.path(.shape_dir,'..','lib','run_manifest.py')); stopifnot(nrow(x$records)==2, all(x$records$input_format=='v2'), grepl('120000',x$run_decisions$selected_run));"
            result = subprocess.run(["Rscript", "-e", expression, str(work / "cohort")], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            cli = subprocess.run(["Rscript", str(ROOT / "src/visualize/spectrum_shape.R"), str(work / "cohort"), "--input-mode", "mixed", "--region", "whole", "--output-dir", str(work / "mixed-output")], capture_output=True, text=True)
            # Fixture spectra may be incomplete, but settings are always audited.
            settings = json.loads((work / "mixed-output/settings.json").read_text())
            self.assertEqual(settings["provenance"], "mixed input; selected v2 runs validated, legacy runs not v2-verified")

            missing = expression.replace("'whole'", "'noncoding'").replace("nrow(x$records)==2,", "nrow(x$matrix)==0, all(!x$records$included), nrow(x$records)==2,")
            result = subprocess.run(["Rscript", "-e", missing, str(work / "cohort")], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            # A single available region artifact leaves the other slot audited.
            import hashlib
            channels = [f"{a}[{v}]{b}" for v in ("C>A","C>G","C>T","T>A","T>C","T>G") for a in "ACGT" for b in "ACGT"]
            partial = run / "partial_ncds.tsv"
            partial.write_text("mutType\tmutNum\ttotalRootNum\n" + "".join(f"{c}\t1\t100\n" for c in channels))
            manifest["artifacts"].append({"kind":"single_ncds", "slot":"org2", "path":partial.name, "sha256":hashlib.sha256(partial.read_bytes()).hexdigest()})
            for i, kind in enumerate(("pdf_ncds_norm","pdf_ncds_sbst","pdf_ncds_ori","pdf_ncds_logratio","pdf_dinuc_ncds_tsv")):
                plot = run / f"partial_plot_{i}.pdf"
                plot.write_bytes(b"%PDF-1.4\nfixture\n%%EOF\n")
                manifest["artifacts"].append({"kind":kind,"slot":"org2","path":plot.name,"sha256":hashlib.sha256(plot.read_bytes()).hexdigest()})
            (run / "metadata/metadata_manifest.json").write_text(json.dumps(manifest))
            partial_check = missing.replace("nrow(x$matrix)==0, all(!x$records$included)", "nrow(x$matrix)==1, sum(!x$records$included)==1")
            result = subprocess.run(["Rscript", "-e", partial_check, str(work / "cohort")], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            duplicate = trio / "bench_20260101_120000_duplicate"
            duplicate.symlink_to(run.resolve(), target_is_directory=True)
            result = subprocess.run(["Rscript", "-e", expression, str(work / "cohort")], cwd=ROOT, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            duplicate.unlink()
            previous = {"organisms": json.loads(json.dumps(manifest["organisms"]))}
            previous["organisms"][0]["accession"] = "GCA_999999999.1"
            (old / "metadata_manifest.json").write_text(json.dumps(previous))
            result = subprocess.run(["Rscript", "-e", expression, str(work / "cohort")], cwd=ROOT, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            (old / "metadata_manifest.json").write_text(json.dumps({"organisms": manifest["organisms"]}))
            manifest["status"] = "pending"
            (run / "metadata/metadata_manifest.json").write_text(json.dumps(manifest))
            result = subprocess.run(["Rscript", "-e", expression, str(work / "cohort")], cwd=ROOT, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            manifest["status"] = "complete"
            # Missing schema_version must still be treated as corrupt v2.
            manifest.pop("schema_version")
            (run / "metadata/metadata_manifest.json").write_text(json.dumps(manifest))
            result = subprocess.run(["Rscript", "-e", expression, str(work / "cohort")], cwd=ROOT, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            manifest["schema_version"] = 2
            (run / "metadata/metadata_manifest.json").write_text(json.dumps(manifest))
            # Corrupt a selected v2 artifact: no fallback to the old legacy run.
            artifact = next(a for a in manifest["artifacts"] if a["kind"] == "single")
            (run / artifact["path"]).write_text("corrupt")
            result = subprocess.run(["Rscript", "-e", expression, str(work / "cohort")], cwd=ROOT, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)

    def test_region_conflict_precedes_input(self):
        result = subprocess.run(["Rscript", str(ROOT / "src/visualize/spectrum_shape.R"), "/missing-input", "--region", "whole", "--tsv-pattern", "*_ncds.tsv"], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Conflicting --region and --tsv-pattern", result.stderr)

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
            colors = read("annotation_colors.csv")
            self.assertEqual(set(colors[0]), {"level", "label", "color", "identifier"})
            self.assertEqual(len({row["color"] for row in colors}), len(colors))
            self.assertEqual(len({row["identifier"] for row in colors}), len(colors))
            for row in colors:
                prefix = {"phylum": "P", "class": "C", "order": "O", "cluster": "K"}.get(row["level"], row["level"])
                self.assertTrue(row["identifier"].startswith(prefix), row)
                if row["level"] == "cluster":
                    self.assertEqual(row["identifier"], "K" + row["label"])
            self.assertEqual({row["level"] for row in colors}, {"phylum", "class", "order", "species", "cluster"})
            for name in ("heatmap.pdf", "heatmap.png", "umap_phylum.png", "umap_cluster.pdf", "spectrum_cluster_1.png"):
                self.assertGreater((work / "output" / name).stat().st_size, 1000)
            # Inspect labels when the local PDF text utility is available.
            import shutil
            if shutil.which("pdftotext"):
                heatmap_text = subprocess.check_output(["pdftotext", "-layout", str(work / "output/heatmap.pdf"), "-"], text=True)
                umap_text = subprocess.check_output(["pdftotext", "-layout", str(work / "output/umap_cluster.pdf"), "-"], text=True)
                for cluster in {row["cluster"] for row in read("clusters.csv")}:
                    self.assertIn("Cluster " + cluster, heatmap_text)
                    self.assertNotIn("K" + cluster + ": " + cluster, heatmap_text)
                    self.assertIn("K" + cluster + ": " + cluster, umap_text)
            # The same counts through legacy ingestion retain numerical results.
            for group in range(3):
                manifest_path = next((work / "input" / str(group)).glob("*/metadata/metadata_manifest.json"))
                manifest = json.loads(manifest_path.read_text())
                legacy = work / "legacy" / str(group) / "bench_20260101"
                (legacy / "statistics/singlenuc").mkdir(parents=True)
                shutil.copytree(manifest_path.parent, legacy / "metadata")
                for artifact in manifest["artifacts"]:
                    if artifact["kind"] == "single":
                        entry = next(e for e in manifest["organisms"] if e["slot"] == artifact["slot"])
                        target = entry["accession"] + "_" + entry["short_name"] + "_test.tsv"
                        shutil.copyfile(manifest_path.parents[1] / artifact["path"], legacy / "statistics/singlenuc" / target)
                (legacy / "metadata/metadata_manifest.json").write_text(json.dumps(manifest))
            legacy_result = subprocess.run(["Rscript", str(ROOT / "src/visualize/spectrum_shape.R"), str(work / "legacy"), "--input-mode", "legacy", "--region", "whole", "--output-dir", str(work / "legacy-output")], capture_output=True, text=True)
            self.assertEqual(legacy_result.returncode, 0, legacy_result.stderr)
            with (work / "legacy-output/log2_matrix.csv").open() as handle:
                legacy_rows = list(csv.reader(handle))[1:]
            with (work / "output/log2_matrix.csv").open() as handle:
                v2_rows = list(csv.reader(handle))[1:]
            self.assertEqual([r[1:] for r in legacy_rows], [r[1:] for r in v2_rows])
            with (work / "legacy-output/clusters.csv").open() as handle:
                legacy_clusters = [r["cluster"] for r in csv.DictReader(handle)]
            self.assertEqual(legacy_clusters, [r["cluster"] for r in read("clusters.csv")])
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
