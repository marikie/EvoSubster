import math
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "statistics"))

import chi2_spectrum_cmp as script

TSV_A = "./fixtures/chi2_cmp_A.tsv"
TSV_B = "./fixtures/chi2_cmp_B.tsv"

# Expected values (computed by hand):
#   For A[C>A]A: r_hat = 300/20000 = 0.015, E_A = E_B = 150
#       contrib = (100-150)^2/150 + (200-150)^2/150 = 100/3
#   For A[C>T]A: same by symmetry -> contrib = 100/3
#   chi2 = 200/3, df = 2, total_mut = 600, V = sqrt((200/3)/600) = 1/3
EXPECTED_CHI2 = 200 / 3
EXPECTED_V = 1 / 3


class TestLoadTsv(unittest.TestCase):
    def test_load(self):
        data = script.load_tsv(TSV_A)
        self.assertIn("A[C>A]A", data)
        self.assertEqual(data["A[C>A]A"], (100.0, 10000.0))
        self.assertEqual(data["A[C>T]A"], (200.0, 10000.0))


class TestComputeChi2(unittest.TestCase):
    def setUp(self):
        self.data_a = script.load_tsv(TSV_A)
        self.data_b = script.load_tsv(TSV_B)

    def test_chi2_stat(self):
        chi2_stat, df, _, _, _ = script.compute_chi2(self.data_a, self.data_b)
        self.assertAlmostEqual(chi2_stat, EXPECTED_CHI2, places=4)
        self.assertEqual(df, 2)

    def test_cramers_v(self):
        _, _, _, cramers_v, _ = script.compute_chi2(self.data_a, self.data_b)
        self.assertAlmostEqual(cramers_v, EXPECTED_V, places=4)

    def test_p_value_significant(self):
        _, _, p_value, _, _ = script.compute_chi2(self.data_a, self.data_b)
        self.assertLess(p_value, 0.05)

    def test_per_type_length(self):
        _, _, _, _, per_type = script.compute_chi2(self.data_a, self.data_b)
        self.assertEqual(len(per_type), 2)

    def test_per_type_contributions_sum(self):
        chi2_stat, _, _, _, per_type = script.compute_chi2(self.data_a, self.data_b)
        total = sum(c for _, c, _ in per_type)
        self.assertAlmostEqual(total, chi2_stat, places=10)

    def test_identical_spectra_chi2_zero(self):
        chi2_stat, _, _, _, _ = script.compute_chi2(self.data_a, self.data_a)
        self.assertAlmostEqual(chi2_stat, 0.0, places=10)

    def test_zero_root_skipped(self):
        data_zero = {"A[C>A]A": (0.0, 0.0)}
        chi2_stat, df, _, _, _ = script.compute_chi2(data_zero, data_zero)
        self.assertEqual(df, 0)
        self.assertTrue(math.isnan(script.compute_chi2(data_zero, data_zero)[2]))


if __name__ == "__main__":
    unittest.main()
