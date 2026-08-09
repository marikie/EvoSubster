import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "select" / "strip_newick_accessions.py"
SPEC = importlib.util.spec_from_file_location("strip_newick_accessions", SCRIPT)
CONVERTER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONVERTER)


class StripNewickAccessionsCliTest(unittest.TestCase):
    def run_converter(self, source: str, *extra_args: str):
        temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)
        work = Path(temp_dir.name)
        input_path = work / "input.nwk"
        output_path = work / "output.nwk"
        input_path.write_text(source, encoding="utf-8")
        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--input",
                str(input_path),
                "--output",
                str(output_path),
                *extra_args,
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        return result, input_path, output_path

    def test_rewrites_only_terminal_accessions_on_leaf_labels(self):
        source = (
            "(([&tip]Chaunax_sp._Z400_GCA_037577475.1[&note=GCA_9.9]:0.1,"
            "A_a[inside]_GCA_000000001.1:0.2)"
            "'Internal_GCA_000000003.3':0.3,"
            "'A''quoted_species_GCF_000001405.40':0.4,Homo_sapiens:0.5);\n"
        )
        expected = (
            "(([&tip]Chaunax_sp._Z400[&note=GCA_9.9]:0.1,"
            "A_a[inside]:0.2)'Internal_GCA_000000003.3':0.3,"
            "'A''quoted_species':0.4,Homo_sapiens:0.5);\n"
        )

        result, _, output_path = self.run_converter(source)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(output_path.read_text(encoding="utf-8"), expected)
        self.assertIn("converted 3 of 4 leaf labels", result.stderr.lower())
        self.assertIn("1 leaf label had no terminal accession", result.stderr.lower())

    def test_leaves_non_versioned_or_species_only_labels_unchanged(self):
        source = "(Homo_sapiens:1,Mus_musculus_GCA_000001635:2);\n"

        result, _, output_path = self.run_converter(source)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(output_path.read_text(encoding="utf-8"), source)
        self.assertIn("converted 0 of 2 leaf labels", result.stderr.lower())
        self.assertIn("2 leaf labels had no terminal accession", result.stderr.lower())

    def test_refuses_to_overwrite_existing_output_without_force(self):
        result, _, output_path = self.run_converter("(A_a_GCA_000000001.1:1,B_b:2);\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        output_path.write_text("keep me\n", encoding="utf-8")

        second = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--input",
                str(output_path.parent / "input.nwk"),
                "--output",
                str(output_path),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertNotEqual(second.returncode, 0)
        self.assertIn("already exists", second.stderr.lower())
        self.assertEqual(output_path.read_text(encoding="utf-8"), "keep me\n")

    def test_force_overwrites_an_existing_output(self):
        result, input_path, output_path = self.run_converter("(A_a_GCA_000000001.1:1,B_b:2);\n")
        self.assertEqual(result.returncode, 0, result.stderr)
        output_path.write_text("stale\n", encoding="utf-8")

        forced = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--input",
                str(input_path),
                "--output",
                str(output_path),
                "--force",
            ],
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(forced.returncode, 0, forced.stderr)
        self.assertEqual(output_path.read_text(encoding="utf-8"), "(A_a:1,B_b:2);\n")

    def test_atomic_writer_refuses_existing_destination_without_overwrite(self):
        temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)
        output_path = Path(temp_dir.name) / "output.nwk"
        output_path.write_text("existing\n", encoding="utf-8")

        with self.assertRaises(FileExistsError):
            CONVERTER.write_text_atomic(output_path, "replacement\n", overwrite=False)

        self.assertEqual(output_path.read_text(encoding="utf-8"), "existing\n")

    def test_refuses_same_input_and_output_path(self):
        temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)
        tree_path = Path(temp_dir.name) / "tree.nwk"
        source = "(A_a_GCA_000000001.1:1,B_b:2);\n"
        tree_path.write_text(source, encoding="utf-8")

        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--input",
                str(tree_path),
                "--output",
                str(tree_path),
                "--force",
            ],
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("must be different", result.stderr.lower())
        self.assertEqual(tree_path.read_text(encoding="utf-8"), source)


class StripNewickAccessionsInterfaceTest(unittest.TestCase):
    def test_scan_leaf_segment_exposes_logical_label_and_source_positions(self):
        text = "A_a[inside]_GCA_000000001.1:0.2"

        scan = CONVERTER.scan_leaf_segment(text, 0)

        self.assertEqual(scan.end, text.index(":"))
        self.assertEqual(scan.logical_label, "A_a_GCA_000000001.1")
        self.assertTrue(scan.has_label)
        self.assertEqual(
            "".join(text[position] for position in scan.source_positions),
            scan.logical_label,
        )

    def test_rewrite_newick_returns_named_counts(self):
        result = CONVERTER.rewrite_newick("(A_a_GCA_000000001.1:1,B_b:2);\n")

        self.assertIsInstance(result, CONVERTER.ConversionResult)
        self.assertEqual(result.text, "(A_a:1,B_b:2);\n")
        self.assertEqual((result.converted, result.leaf_count), (1, 2))

    def test_convert_file_owns_path_validation(self):
        temp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(temp_dir.cleanup)
        tree_path = Path(temp_dir.name) / "tree.nwk"
        tree_path.write_text("(A_a_GCA_000000001.1:1,B_b:2);\n", encoding="utf-8")

        with self.assertRaisesRegex(
            CONVERTER.CliUsageError,
            "--input and --output must be different paths",
        ):
            CONVERTER.convert_file(tree_path, tree_path, overwrite=True)


if __name__ == "__main__":
    unittest.main()
