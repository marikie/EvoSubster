import os
import sys
import unittest
import subprocess

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "count"))

import single_sbst_2TSVs_errprob as script


class TestTriUvMuts2TSVs(unittest.TestCase):
    def test_main(self):
        print("test_main")

        # inputs
        self.in1_maf = open("./fixtures/test_triUvMuts2TSVs_errprb_in1.maf")
        self.org2outPath1 = "./fixtures/result_triUvMuts2TSVs_errprb_out1_2.tsv"
        self.org3outPath1 = "./fixtures/result_triUvMuts2TSVs_errprb_out1_3.tsv"

        # expected outputs
        self.out1_2 = "./fixtures/test_triUvMuts2TSVs_errprb_out1_2.tsv"
        self.out1_3 = "./fixtures/test_triUvMuts2TSVs_errprb_out1_3.tsv"

        ### test1 ###
        print("test1")
        # Run the main function with test data
        script.main(self.in1_maf, self.org2outPath1, self.org3outPath1)
        # Check if the contents of the two files (result and expected) are the same using diff
        result = subprocess.run(
            ["diff", "-b", self.org2outPath1, self.out1_2], capture_output=True
        )
        self.assertEqual(result.returncode, 0, "The files are different")
        result = subprocess.run(
            ["diff", "-b", self.org3outPath1, self.out1_3], capture_output=True
        )
        self.assertEqual(result.returncode, 0, "The files are different")

    def test_bothEdgeBasesSame(self):
        print("test_bothEdgeBasesSame")
        self.assertTrue(script.bothEdgeBasesSame("AAA", "ATA", "ACA"))
        self.assertTrue(script.bothEdgeBasesSame("AAA", "AAA", "AAA"))
        self.assertFalse(script.bothEdgeBasesSame("TAA", "ATA", "AGA"))
        self.assertFalse(script.bothEdgeBasesSame("AAA", "AAA", "AAC"))

    def test_rev(self):
        print("test_rev")
        self.assertEqual(script.rev("AAA"), "TTT")
        self.assertEqual(script.rev("ATA"), "TAT")
        self.assertEqual(script.rev("ACA"), "TGT")
        self.assertEqual(script.rev("ACC"), "GGT")

    def test_mutType(self):
        print("test_mutType")
        self.assertEqual(script.mutType("AAA", "A", "G"), "TTTC")
        self.assertEqual(script.mutType("ACG", "C", "T"), "ACGT")

    def test_add2totalNum(self):
        print("test_add2totalNum")
        mutDict = script.initialize_mut_dict()
        exp_mutDict = script.initialize_mut_dict()
        exp_mutDict["CTGA"] = {"mutNum": 0, "totalRootNum": 1}
        exp_mutDict["CTGC"] = {"mutNum": 0, "totalRootNum": 1}
        exp_mutDict["CTGG"] = {"mutNum": 0, "totalRootNum": 1}
        script.add2totalNum(mutDict, "CTG")
        self.assertEqual(mutDict, exp_mutDict)

        mutDict = script.initialize_mut_dict()
        exp_mutDict = script.initialize_mut_dict()
        exp_mutDict["CTGA"] = {"mutNum": 0, "totalRootNum": 1}
        exp_mutDict["CTGC"] = {"mutNum": 0, "totalRootNum": 1}
        exp_mutDict["CTGG"] = {"mutNum": 0, "totalRootNum": 1}
        script.add2totalNum(mutDict, "CAG")
        self.assertEqual(mutDict, exp_mutDict)

    def test_add2MutDict(self):
        print("test_add2MutDict")
        mutDict2 = script.initialize_mut_dict()
        mutDict3 = script.initialize_mut_dict()
        script.add2MutDict("CAG", "CGG", "CAG", mutDict2, mutDict3)
        exp_mutDict2 = script.initialize_mut_dict()
        exp_mutDict2["CTGA"] = {"mutNum": 0, "totalRootNum": 1}
        exp_mutDict2["CTGC"] = {"mutNum": 1, "totalRootNum": 1}
        exp_mutDict2["CTGG"] = {"mutNum": 0, "totalRootNum": 1}
        exp_mutDict3 = script.initialize_mut_dict()
        exp_mutDict3["CTGA"] = {"mutNum": 0, "totalRootNum": 1}
        exp_mutDict3["CTGC"] = {"mutNum": 0, "totalRootNum": 1}
        exp_mutDict3["CTGG"] = {"mutNum": 0, "totalRootNum": 1}

        # Set maxDiff to None to see the full diff
        self.maxDiff = None

        self.assertEqual(mutDict2, exp_mutDict2)
        self.assertEqual(mutDict3, exp_mutDict3)


class TestPathOperations(unittest.TestCase):
    def setUp(self):
        # Test cases that will be used across multiple tests
        self.test_paths = [
            "/home/user/data/ulvPro_ulvMut_ulvCom_aligned.maf",
            "data/ulvPro_ulvMut_ulvCom_aligned.maf",
            "./ulvPro_ulvMut_ulvCom_aligned.maf",
            "ulvPro_ulvMut_ulvCom_aligned.maf",
        ]

    def test_filename_extraction(self):
        expected_filename = "ulvPro_ulvMut_ulvCom_aligned"
        for path in self.test_paths:
            with self.subTest(path=path):
                filename = os.path.splitext(os.path.basename(path))[0]
                self.assertEqual(filename, expected_filename)

    def test_directory_path(self):
        expected_dirs = ["/home/user/data", "data", ".", ""]
        for path, expected_dir in zip(self.test_paths, expected_dirs):
            with self.subTest(path=path):
                path_before_filename = os.path.dirname(path)
                self.assertEqual(path_before_filename, expected_dir)

    def test_filename_parts(self):
        expected_parts = ["ulvPro", "ulvMut", "ulvCom", "aligned"]
        for path in self.test_paths:
            with self.subTest(path=path):
                filename = os.path.splitext(os.path.basename(path))[0]
                filename_parts = filename.split("_")
                self.assertEqual(filename_parts, expected_parts)

    def test_output_path_construction(self):
        expected_outputs = [
            (
                "/home/user/data/ulvMut_aligned.tsv",
                "/home/user/data/ulvCom_aligned.tsv",
            ),
            ("data/ulvMut_aligned.tsv", "data/ulvCom_aligned.tsv"),
            ("./ulvMut_aligned.tsv", "./ulvCom_aligned.tsv"),
            ("ulvMut_aligned.tsv", "ulvCom_aligned.tsv"),
        ]

        for path, expected in zip(self.test_paths, expected_outputs):
            with self.subTest(path=path):
                # Get filename without extension
                filename = os.path.splitext(os.path.basename(path))[0]
                path_before_filename = os.path.dirname(path)

                # Split filename
                filename_parts = filename.split("_")
                org2 = filename_parts[1]
                org3 = filename_parts[2]
                rest = "_".join(filename_parts[3:])

                # Construct output paths
                outFile2 = f"{org2}_{rest}.tsv"
                outFile3 = f"{org3}_{rest}.tsv"
                outputFilePath2 = os.path.join(path_before_filename, outFile2)
                outputFilePath3 = os.path.join(path_before_filename, outFile3)

                self.assertEqual((outputFilePath2, outputFilePath3), expected)

    def test_invalid_filename(self):
        invalid_paths = [
            "/home/user/data/invalid.maf",
            "data/ulvPro_ulvMut.maf",
            "not_enough_parts.maf",
        ]

        for path in invalid_paths:
            with self.subTest(path=path):
                filename = os.path.splitext(os.path.basename(path))[0]
                filename_parts = filename.split("_")

                # First check if we have enough parts
                if len(filename_parts) < 3:
                    with self.assertRaises(ValueError):
                        self.validate_filename_parts(filename_parts)

    def validate_filename_parts(self, parts):
        if len(parts) < 3:
            raise ValueError("The file name should be org1_org2_org3_*.maf")
        return parts[0], parts[1], parts[2]


if __name__ == "__main__":
    unittest.main()
