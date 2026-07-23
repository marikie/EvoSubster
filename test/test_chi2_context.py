import math
import os
import re
import shutil
import subprocess
import tempfile
import unittest


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT_PATH = os.path.join(ROOT_DIR, "src", "statistics", "chi2_context.R")
RSCRIPT = shutil.which("Rscript")


@unittest.skipUnless(RSCRIPT, "Rscript is required for chi2_context.R tests")
class Chi2ContextRTest(unittest.TestCase):
    def run_script(self, rows, label=""):
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = os.path.join(tmp_dir, "contexts.tsv")
            with open(input_path, "w", encoding="utf-8") as handle:
                handle.write("mutType\tmutNum\ttotalRootNum\n")
                for mut_type, mut_num, root_num in rows:
                    handle.write(f"{mut_type}\t{mut_num}\t{root_num}\n")

            command = [RSCRIPT, SCRIPT_PATH, input_path]
            if label:
                command.append(label)
            return subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
            )

    def parse_result(self, stdout, mutation_type):
        pattern = re.compile(
            rf"^{re.escape(mutation_type)}\s+"
            r"(?P<method>\w+)\s+"
            r"(?P<statistic>\S+)\s+"
            r"(?P<df>\S+)\s+"
            r"(?P<p_value>\S+)",
            re.MULTILINE,
        )
        match = pattern.search(stdout)
        self.assertIsNotNone(match, stdout)
        return match.groupdict()

    def test_chisq_result_matches_contingency_table(self):
        rows = [
            ("A[C>A]A", 100, 10000),
            ("C[C>A]A", 50, 10000),
            ("G[C>A]A", 100, 10000),
            ("T[C>A]A", 50, 10000),
        ]

        result = self.run_script(rows, label="synthetic")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("synthetic", result.stdout)
        parsed = self.parse_result(result.stdout, "C>A")
        self.assertEqual(parsed["method"], "chisq")
        self.assertEqual(int(parsed["df"]), 3)

        expected_mut = 75.0
        expected_non = 9925.0
        expected_chi2 = sum(
            ((mut - expected_mut) ** 2 / expected_mut)
            + (((root - mut) - expected_non) ** 2 / expected_non)
            for _, mut, root in rows
        )
        self.assertAlmostEqual(
            float(parsed["statistic"]),
            expected_chi2,
            places=4,
        )
        self.assertLess(float(parsed["p_value"]), 0.05)

    def test_uniform_context_rates_give_zero_chi_square(self):
        rows = [
            ("A[C>A]A", 75, 10000),
            ("C[C>A]A", 75, 10000),
            ("G[C>A]A", 75, 10000),
            ("T[C>A]A", 75, 10000),
        ]

        result = self.run_script(rows)

        self.assertEqual(result.returncode, 0, result.stderr)
        parsed = self.parse_result(result.stdout, "C>A")
        self.assertEqual(parsed["method"], "chisq")
        self.assertTrue(math.isclose(float(parsed["statistic"]), 0.0, abs_tol=1e-10))

    def test_rejects_input_without_positive_root_counts(self):
        result = self.run_script([("A[C>A]A", 0, 0)])

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("No rows with totalRootNum > 0", result.stderr)


if __name__ == "__main__":
    unittest.main()
