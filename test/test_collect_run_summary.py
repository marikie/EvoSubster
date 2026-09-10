import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "report"))

import collect_run_summary as script
import report_contract


SLOTS = (
    ("org1", "Out1", "GCF_1", "Pan alpha"),
    ("org2", "In2", "GCA_2", "Homo beta"),
    ("org3", "In3", "GCA_3", "Homo gamma"),
)


def write_file(path, content=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def create_run(dataset_dir, layout, date="20260101", stale_metadata=False):
    run_dir = dataset_dir / date
    metadata_dir = run_dir / "metadata"
    organisms = []

    for slot, short_name, accession, organism_name in SLOTS:
        metadata_path = metadata_dir / f"{short_name}_{accession}.json"
        write_file(metadata_path, "{}")
        write_file(
            metadata_dir / f"taxonomy_{metadata_path.name}",
            json.dumps(
                {
                    "reports": [
                        {
                            "taxonomy": {
                                "classification": {
                                    "genus": {"name": organism_name.split()[0]},
                                    "species": {"name": organism_name},
                                }
                            }
                        }
                    ]
                }
            ),
        )
        manifest_metadata = (
            Path("/removed/results") / dataset_dir.name / date / "metadata" / metadata_path.name
            if stale_metadata
            else metadata_path
        )
        organisms.append(
            {
                "slot": slot,
                "short_name": short_name,
                "accession": accession,
                "raw_organism_name": organism_name,
                "metadata_json": str(manifest_metadata),
            }
        )

    write_file(metadata_dir / "metadata_manifest.json", json.dumps({"organisms": organisms}))

    artifact_root = run_dir if layout == "legacy" else run_dir / "statistics" / "misc"
    write_file(artifact_root / f"Out1_gcContent_{date}.out", "Total GC content: 40%\n")
    write_file(artifact_root / f"In2_gcContent_{date}.out", "Total GC content: 41%\n")
    write_file(artifact_root / f"In3_gcContent_{date}.out", "Total GC content: 42%\n")
    write_file(artifact_root / f"sbstRatio_{date}.out", "org2: 1.2\norg3: 1.3\n")

    train_root = run_dir if layout == "legacy" else run_dir / "intermediateFiles"
    for pair, identity in (("Out12In2", 90), ("Out12In3", 91), ("In22In3", 97)):
        write_file(
            train_root / f"{pair}_{date}.train",
            f"# substitution percent identity: {identity}\n",
        )

    for _, short_name, accession, _ in SLOTS[1:]:
        if layout == "legacy":
            tsv_dir = run_dir
            pdf_dir = run_dir
        else:
            tsv_dir = run_dir / "statistics" / short_name / "singlenuc"
            pdf_dir = run_dir / "figs" / short_name / "singlenuc"
        write_file(tsv_dir / f"{accession}_{short_name}_{date}.tsv", "context\tcount\n")
        write_file(
            pdf_dir / ("ratio" if layout == "current" else "") / f"{accession}_{short_name}_{date}_norm.pdf"
        )
        write_file(
            pdf_dir / ("log-ratio" if layout == "current" else "") / f"{accession}_{short_name}_{date}_logRatio_mean.pdf"
        )
        write_file(
            (run_dir / "figs" / short_name / "dinuc" if layout == "current" else pdf_dir)
            / f"{accession}_{short_name}_{date}_dinuc.tsv.pdf"
        )

    return run_dir


class BuildSummaryCompatibilityTest(unittest.TestCase):
    def test_reads_current_layout_and_rebases_stale_metadata_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            dataset_dir = Path(tmp) / "TrioA"
            run_dir = create_run(dataset_dir, "current", stale_metadata=True)

            summary, _ = script.build_summary(dataset_dir, 80.0)

            self.assertEqual(len(summary["datasets"]), 1)
            entry = summary["datasets"][0]
            self.assertEqual(entry["species2"]["metadata_json"], str(run_dir / "metadata" / "In2_GCA_2.json"))
            self.assertIn("statistics/In2/singlenuc", entry["species2"]["tsv_file"])
            self.assertEqual(summary["issues"], [])

    def test_prefers_local_metadata_even_when_original_path_still_exists(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dataset_dir = root / "TrioA"
            run_dir = create_run(dataset_dir, "current")
            external_metadata = root / "original" / "In2_GCA_2.json"
            write_file(external_metadata, "{}")
            manifest_path = run_dir / "metadata" / "metadata_manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["organisms"][1]["metadata_json"] = str(external_metadata)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            summary, _ = script.build_summary(dataset_dir, 80.0)

            self.assertEqual(
                summary["datasets"][0]["species2"]["metadata_json"],
                str(run_dir / "metadata" / "In2_GCA_2.json"),
            )

    def test_reads_legacy_run_root_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            dataset_dir = Path(tmp) / "TrioLegacy"
            run_dir = create_run(dataset_dir, "legacy")

            summary, _ = script.build_summary(dataset_dir, 80.0)

            entry = summary["datasets"][0]
            self.assertEqual(Path(entry["species2"]["tsv_file"]).parent, run_dir)
            self.assertEqual(entry["idt_23"], "97 %")
            self.assertEqual(summary["issues"], [])

    def test_prefers_current_artifacts_when_legacy_files_also_exist(self):
        with tempfile.TemporaryDirectory() as tmp:
            dataset_dir = Path(tmp) / "TrioMixed"
            run_dir = create_run(dataset_dir, "current")
            write_file(run_dir / "GCA_2_In2_20260101.tsv", "legacy\n")

            summary, _ = script.build_summary(dataset_dir, 80.0)

            self.assertIn("statistics/In2/singlenuc", summary["datasets"][0]["species2"]["tsv_file"])

    def test_pdf_search_falls_back_after_excluding_current_ncds_only_matches(self):
        with tempfile.TemporaryDirectory() as tmp:
            dataset_dir = Path(tmp) / "TrioMixed"
            run_dir = create_run(dataset_dir, "current")
            for short_name, accession in (("In2", "GCA_2"), ("In3", "GCA_3")):
                (run_dir / "figs" / short_name / "singlenuc" / "ratio" / f"{accession}_{short_name}_20260101_norm.pdf").unlink()
                (run_dir / "figs" / short_name / "singlenuc" / "log-ratio" / f"{accession}_{short_name}_20260101_logRatio_mean.pdf").unlink()
                write_file(run_dir / "figs" / short_name / "singlenuc" / "ratio" / f"{accession}_{short_name}_20260101_ncds_norm.pdf")
                write_file(run_dir / "figs" / short_name / "singlenuc" / "log-ratio" / f"{accession}_{short_name}_20260101_ncds_logRatio_mean.pdf")
                write_file(run_dir / f"{accession}_{short_name}_20260101_norm.pdf")
                write_file(run_dir / f"{accession}_{short_name}_20260101_logRatio_mean.pdf")

            summary, _ = script.build_summary(dataset_dir, 80.0)

            self.assertEqual(
                Path(summary["datasets"][0]["species2"]["pdfs"]["norm"][0]).parent,
                run_dir,
            )
            self.assertEqual(summary["issues"], [])

    def test_lineage_root_ignores_non_dataset_directories(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "lineage"
            create_run(root / "TrioA", "current")
            (root / "pca" / "20260101").mkdir(parents=True)
            (root / "tmp").mkdir()

            summary, _ = script.build_summary(root, 80.0)

            self.assertEqual([entry["dataset"] for entry in summary["datasets"]], ["TrioA"])
            self.assertEqual(summary["issues"], [])

    def test_reports_missing_manifest_in_newest_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "lineage"
            dataset_dir = root / "TrioA"
            create_run(dataset_dir, "current")
            (dataset_dir / "20260102").mkdir()

            summary, _ = script.build_summary(root, 80.0)

            self.assertEqual(summary["datasets"], [])
            self.assertEqual(len(summary["issues"]), 1)
            self.assertIn("Metadata manifest not found", summary["issues"][0])

    def test_reports_missing_manifest_in_first_incomplete_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "lineage"
            (root / "Out1_In2_In3" / "20260101").mkdir(parents=True)
            (root / "pca" / "20260101").mkdir(parents=True)

            summary, _ = script.build_summary(root, 80.0)

            self.assertEqual(summary["datasets"], [])
            self.assertEqual(len(summary["issues"]), 1)
            self.assertTrue(summary["issues"][0].startswith("Out1_In2_In3:"))
            self.assertIn("Metadata manifest not found", summary["issues"][0])


class ReportContractTest(unittest.TestCase):
    def test_derives_default_and_filtered_summary_paths(self):
        root = Path("/tmp/results/fungi")
        summary, filtered = report_contract.summary_paths(root)
        self.assertEqual(summary, root / "fungi_summary.json")
        self.assertEqual(filtered, root / "fungi_summary_filtered.json")

    def test_report_extension_matches_output_format(self):
        summary = Path("/tmp/fungi_summary.json")
        self.assertEqual(
            report_contract.report_output_path(summary, "word_document"),
            Path("/tmp/fungi_summary.docx"),
        )
        self.assertEqual(
            report_contract.report_output_path(summary, "html_document"),
            Path("/tmp/fungi_summary.html"),
        )
        self.assertEqual(
            report_contract.report_output_path(summary, "pdf_document"),
            Path("/tmp/fungi_summary.pdf"),
        )

    def test_rejects_unsupported_output_format(self):
        with self.assertRaisesRegex(ValueError, "unsupported output format"):
            report_contract.report_output_path(
                Path("/tmp/fungi_summary.json"), "github_document"
            )

    def test_validates_required_dataset_fields(self):
        dataset = {"dataset": "TrioA"}
        names, errors = report_contract.validate_summary({"datasets": [dataset]})
        self.assertEqual(names, ["TrioA"])
        self.assertEqual(len(errors), 1)
        self.assertIn("species1", errors[0])
        self.assertIn("idt_23", errors[0])

    def test_rejects_invalid_species_and_pdf_types(self):
        dataset = {
            "dataset": "TrioA",
            "species1": {"metadata": "not-an-object"},
            "species2": "not-an-object",
            "species3": {
                "taxonomy": "not-an-object",
                "gc_content": "not-an-array",
                "pdfs": {"norm": "not-an-array"},
            },
            "idt_12": "90 %",
            "idt_13": 91,
            "idt_23": "97 %",
        }

        _, errors = report_contract.validate_summary({"datasets": [dataset]})

        self.assertTrue(any("species2 must be an object" in error for error in errors))
        self.assertTrue(any("metadata must be an object" in error for error in errors))
        self.assertTrue(any("taxonomy must be an object" in error for error in errors))
        self.assertTrue(any("gc_content must be an array" in error for error in errors))
        self.assertTrue(any("pdfs.norm must be an array" in error for error in errors))
        self.assertTrue(any("idt_13 must be a string" in error for error in errors))


if __name__ == "__main__":
    unittest.main()