import os
import sys
import unittest

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


if __name__ == "__main__":
    unittest.main()
