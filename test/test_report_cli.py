import json
import shutil
import subprocess
import tempfile
import unittest
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
RUN_REPORT = ROOT / "src" / "report" / "run_report.sh"
RENDER_REPORT = ROOT / "src" / "report" / "render_report.sh"


class ReportCliTest(unittest.TestCase):
    @unittest.skipUnless(shutil.which("Rscript") and shutil.which("pandoc"),
                         "DOCX rendering requires Rscript and pandoc")
    def test_word_report_contains_coding_percentages_and_missing_placeholder(self):
        dependencies = self.run_command(
            "Rscript", "--vanilla", "-e",
            'p <- c("rmarkdown", "jsonlite", "dplyr", "purrr", "tidyr", '
            '"stringr", "glue", "knitr"); '
            'quit(status=if(all(vapply(p, requireNamespace, logical(1), quietly=TRUE))) 0 else 1)',
        )
        if dependencies.returncode:
            self.skipTest("R report packages unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            summary = Path(tmp) / "summary.json"
            output = Path(tmp) / "summary.docx"
            summary.write_text(json.dumps({"datasets": [{
                "dataset": "CodingTest", "latest_run": "20260101",
                "species1": {"coding_region_percentage": ["Coding region percentage: 25.00%"]},
                "species2": {"coding_region_percentage": ["Coding region percentage: 12.34%"]},
                "species3": {}, "idt_12": "90%", "idt_13": "91%", "idt_23": "97%",
            }]}))
            result = self.run_command("bash", str(RENDER_REPORT), "-j", str(summary),
                                      "-o", str(output))
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with zipfile.ZipFile(output) as archive:
                document = ET.fromstring(archive.read("word/document.xml"))
            ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
            table = document.find(".//w:tbl", ns)
            rows = [["".join(cell.itertext()) for cell in row.findall("w:tc", ns)]
                    for row in table.findall("w:tr", ns)]
            self.assertIn("Coding Region (%)", rows[0])
            self.assertEqual([row[-1] for row in rows[1:]], ["25.00%", "12.34%", "—"])
            self.assertIn("including N", "".join(document.itertext()))

    def run_command(self, *args):
        return subprocess.run(
            args,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )

    def test_run_report_help_describes_supported_input_roots(self):
        result = self.run_command("bash", str(RUN_REPORT), "--help")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Lineage root", result.stdout)
        self.assertIn("single trio root", result.stdout)

    def test_run_report_rejects_extra_input_directory(self):
        result = self.run_command("bash", str(RUN_REPORT), ".", ".")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("exactly one input_dir", result.stderr)

    def test_run_report_rejects_empty_equals_value(self):
        result = self.run_command("bash", str(RUN_REPORT), ".", "--json=")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--json requires a value", result.stderr)

    def test_run_report_rejects_unsupported_format_before_collection(self):
        result = self.run_command(
            "bash", str(RUN_REPORT), ".", "--format", "github_document"
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("unsupported output format", result.stderr)
        self.assertNotIn("Step 1", result.stdout)

    def test_run_report_allows_stdout_json_in_collect_only_mode(self):
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "20260101").mkdir()
            result = self.run_command(
                "bash",
                str(RUN_REPORT),
                tmp,
                "--json",
                "-",
                "--collect-only",
            )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No v2 runs found", result.stderr)

    def test_render_report_rejects_missing_option_value(self):
        result = self.run_command("bash", str(RENDER_REPORT), "--json")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires a value", result.stderr)

    def test_render_report_does_not_consume_an_option_as_a_value(self):
        result = self.run_command(
            "bash", str(RENDER_REPORT), "--format", "--help"
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--format requires a value", result.stderr)

    def test_render_report_rejects_invalid_summary_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary_path = Path(tmp) / "invalid.json"
            summary_path.write_text(json.dumps({"datasets": []}), encoding="utf-8")

            result = self.run_command(
                "bash", str(RENDER_REPORT), "--json", str(summary_path)
            )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("summary JSON contains no datasets", result.stderr)


if __name__ == "__main__":
    unittest.main()
