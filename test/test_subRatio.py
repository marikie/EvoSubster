import os
import sys
import unittest
import tempfile
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "metrics"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "count"))

import subRatio as script


class TestSubRatio(unittest.TestCase):
    def test_count(self):
        print("test_count")
        # input
        self.in1 = "./fixtures/test_subRatio_in1.maf"
        result = script.count(self.in1)
        self.assertEqual(result, (5, 11, 4, 12))

    def test_ambiguous_bases_do_not_affect_numerator_or_denominator(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).resolve().parents[1]/'tmp/codex') as tmp:
            path=Path(tmp)/'ambiguous.maf'
            path.write_text('a\ns A 0 5 + 5 ACNNN\ns B 0 5 + 5 ATANN\ns C 0 5 + 5 ACNNN\n')
            self.assertEqual(script.count(path),(1,1,0,2))

    def test_zero_canonical_columns_cannot_produce_percentage(self):
        with self.assertRaises(ValueError):
            script.printPercentage(0,0,0,0)


if __name__ == "__main__":
    unittest.main()
