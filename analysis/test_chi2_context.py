import math
import unittest

import chi2_context as script

TSV_CTX = "./test/chi2_ctx.tsv"

# Expected values for the C>A group (computed by hand):
#   r_hat = (100+50+100+50) / (4*10000) = 300/40000 = 0.0075
#   E for each context = 75
#   contrib per context = 25^2/75 = 625/75
#   chi2 = 4 * 625/75 = 100/3, df = 3
#   eta2 = (100/3) / (100/3 + 300) = 0.1
EXPECTED_CHI2 = 100 / 3
EXPECTED_ETA2 = 0.1


class TestCentralMut(unittest.TestCase):
    def test_parse_standard(self):
        self.assertEqual(script.central_mut("A[C>A]A"), "C>A")
        self.assertEqual(script.central_mut("T[T>G]C"), "T>G")
        self.assertEqual(script.central_mut("G[C>T]T"), "C>T")


class TestGroupByCentral(unittest.TestCase):
    def test_grouping(self):
        data = script.load_tsv(TSV_CTX)
        groups = script.group_by_central(data)
        self.assertIn("C>A", groups)
        self.assertEqual(len(groups["C>A"]), 4)

    def test_only_one_central_type(self):
        data = script.load_tsv(TSV_CTX)
        groups = script.group_by_central(data)
        self.assertEqual(list(groups.keys()), ["C>A"])


class TestComputeChi2Context(unittest.TestCase):
    def setUp(self):
        data = script.load_tsv(TSV_CTX)
        self.groups = script.group_by_central(data)

    def test_chi2_stat(self):
        chi2_stat, df, _, _ = script.compute_chi2_context(self.groups["C>A"])
        self.assertAlmostEqual(chi2_stat, EXPECTED_CHI2, places=4)
        self.assertEqual(df, 3)

    def test_eta2(self):
        _, _, _, eta2 = script.compute_chi2_context(self.groups["C>A"])
        self.assertAlmostEqual(eta2, EXPECTED_ETA2, places=4)

    def test_p_value_significant(self):
        _, _, p_value, _ = script.compute_chi2_context(self.groups["C>A"])
        self.assertLess(p_value, 0.05)

    def test_uniform_rates_give_zero_chi2(self):
        uniform = {
            "A[C>A]A": (75.0, 10000.0),
            "C[C>A]A": (75.0, 10000.0),
            "G[C>A]A": (75.0, 10000.0),
            "T[C>A]A": (75.0, 10000.0),
        }
        chi2_stat, _, _, _ = script.compute_chi2_context(uniform)
        self.assertAlmostEqual(chi2_stat, 0.0, places=10)

    def test_empty_contexts_returns_nan(self):
        empty = {"A[C>A]A": (0.0, 0.0)}
        chi2_stat, df, p_value, eta2 = script.compute_chi2_context(empty)
        self.assertTrue(math.isnan(chi2_stat))
        self.assertEqual(df, 0)


if __name__ == "__main__":
    unittest.main()
