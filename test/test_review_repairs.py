import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('links', ROOT / 'src/report/link_tsvs.py')
links = importlib.util.module_from_spec(spec)
spec.loader.exec_module(links)


class LinkSafetyTests(unittest.TestCase):
    def test_self_directory_preserves_source(self):
        with tempfile.TemporaryDirectory(dir=ROOT / 'tmp/codex') as tmp:
            root = Path(tmp)
            source = root / 'x.tsv'
            source.write_text('source')
            with self.assertRaises(ValueError):
                links.process_links(root, root, '*.tsv')
            self.assertEqual(source.read_text(), 'source')
            self.assertFalse(source.is_symlink())

    def test_preflight_conflict_does_not_create_output(self):
        with tempfile.TemporaryDirectory(dir=ROOT / 'tmp/codex') as tmp:
            root = Path(tmp)
            for name in ('a', 'b'):
                (root / name).mkdir()
                (root / name / 'x.tsv').write_text(name)
            output = root / 'output'
            with self.assertRaises(ValueError):
                links.process_links(root, output, '*.tsv')
            self.assertFalse(output.exists())

    def test_nested_output_and_idempotence(self):
        with tempfile.TemporaryDirectory(dir=ROOT / 'tmp/codex') as tmp:
            root = Path(tmp)
            source = root / 'x.tsv'
            source.write_text('source')
            output = root / 'output'
            links.process_links(root, output, '*.tsv')
            links.process_links(root, output, '*.tsv')
            self.assertEqual((output / 'x.tsv').resolve(), source)

    def test_regular_destination_is_preserved(self):
        with tempfile.TemporaryDirectory(dir=ROOT / 'tmp/codex') as tmp:
            root = Path(tmp)
            (root / 'x.tsv').write_text('source')
            output = root / 'output'
            output.mkdir()
            (output / 'x.tsv').write_text('valuable')
            with self.assertRaises(ValueError):
                links.process_links(root, output, '*.tsv')
            self.assertEqual((output / 'x.tsv').read_text(), 'valuable')
