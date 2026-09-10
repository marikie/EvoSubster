import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
RUN_REPORT = ROOT / "src" / "report" / "run_report.sh"
RENDER_REPORT = ROOT / "src" / "report" / "render_report.sh"


class ReportCliTest(unittest.TestCase):
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

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('"datasets": []', result.stdout)
        self.assertIn("collect-only mode: done", result.stdout)

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
