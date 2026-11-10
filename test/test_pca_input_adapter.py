import json
from pathlib import Path
import sys
import tempfile
import unittest

from test_collect_run_summary import create_run, script

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src/lib'))
import pca_input_adapter as adapter


def create_legacy(root, name='Out1_In2_In3', date='20251027', nested=True):
    run = create_run(root / name, date=date)
    manifest = script.runs.load(run)
    for entry in manifest['organisms']:
        short, accession = entry['short_name'], entry['accession']
        for field, filename in (
            ('metadata_json', f'{short}_{accession}.json'),
            ('taxonomy_json', f'taxonomy_{short}_{accession}.json'),
        ):
            (run / entry[field]).rename(run / 'metadata' / filename)
            entry[field] = '/stale/absolute/path/' + filename
    for artifact in manifest['artifacts']:
        if artifact['kind'] == 'train':
            ref, query = artifact['slot'].split('-')
            names = {entry['slot']: entry['short_name'] for entry in manifest['organisms']}
            target = run / ('intermediateFiles' if nested else '.') / f'{names[ref]}2{names[query]}_{date}.train'
            target.parent.mkdir(exist_ok=True)
            (run / artifact['path']).rename(target)
        if artifact['kind'] == 'single':
            entry = next(e for e in manifest['organisms'] if e['slot'] == artifact['slot'])
            folder = run / 'statistics' / entry['short_name'] / 'singlenuc' if nested else run
            folder.mkdir(parents=True, exist_ok=True)
            filename = f"{entry['accession']}_{entry['short_name']}_{date}"
            baseline_name = 'Oikalb2' if entry['slot'] == 'org2' else 'Oikvan3'
            source = next(p for p in (ROOT / 'test/fixtures/pipeline_counts/20261031101732' / baseline_name / 'singlenuc').glob('*.tsv') if '_ncds' not in p.name)
            (folder / (filename + '.tsv')).write_bytes(source.read_bytes())
            (folder / (filename + '_ncds.tsv')).write_bytes(source.read_bytes())
            for suffix in ('_maflinked.tsv', '_dinuc.tsv'):
                (folder / (filename + suffix)).write_text('invalid derived data')
    legacy = dict(date=date, organisms=manifest['organisms'])
    (run / adapter.runs.MANIFEST).write_text(json.dumps(legacy))
    return run


class PCAInputAdapterTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / 'tmp/codex', prefix='pca-adapter-')
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_flat_and_nested_normalization(self):
        for nested in (False, True):
            root = self.root / str(nested)
            run = create_legacy(root, nested=nested)
            manifest = adapter.read(run)
            self.assertEqual(manifest['trio_key'], 'Out1_In2_In3')
            self.assertEqual(len([a for a in manifest['artifacts'] if a['kind'] == 'single']), 2)
            self.assertEqual(len([a for a in manifest['artifacts'] if a['kind'] == 'single_ncds']), 2)
            self.assertEqual([a['identity_percent'] for a in manifest['artifacts'] if a['kind'] == 'train'], [90, 91, 97])
            for entry in manifest['organisms']:
                self.assertTrue((run / entry['metadata_json']).is_file())
                self.assertTrue((run / entry['taxonomy_json']).is_file())
            self.assertEqual(adapter.discover(root), [str(run)])

    def test_latest_date_symlinks_exclude_analysis_and_bench(self):
        create_legacy(self.root)
        latest = create_legacy(self.root, date='20261031101732')
        create_run(latest.parent / 'bench_run')
        create_run(self.root / 'pca_mixed')
        create_run(self.root / 'spectrum_shape_bench_noncoding')
        create_run(self.root / 'spectrum_shape_bench_whole')
        create_legacy(self.root / 'archive', name='Hidden1_Hidden2_Hidden3')
        collection = self.root / 'collection'
        collection.mkdir()
        (collection / latest.parent.name).symlink_to(latest.parent)
        self.assertEqual(adapter.discover(collection), [str(latest)])
        self.assertEqual(adapter.discover(self.root), [str(latest)])

    def test_ambiguous_canonical_tsv_fails(self):
        run = create_legacy(self.root)
        entry = adapter.read(run)
        artifact = next(a for a in entry['artifacts'] if a['kind'] == 'single')
        (run / Path(artifact['path']).name).write_bytes((run / artifact['path']).read_bytes())
        with self.assertRaisesRegex(ValueError, 'Ambiguous canonical TSV'):
            adapter.read(run)

    def test_missing_inputs_are_visible_and_unsupported_schema_fails(self):
        run = create_legacy(self.root)
        manifest = adapter.read(run)
        artifact = next(a for a in manifest['artifacts'] if a['kind'] == 'single')
        (run / artifact['path']).unlink()
        manifest = adapter.read(run)
        self.assertEqual(len([a for a in manifest['artifacts'] if a['kind'] == 'single']), 1)
        raw = json.loads((run / adapter.runs.MANIFEST).read_text())
        raw['schema_version'] = 1
        (run / adapter.runs.MANIFEST).write_text(json.dumps(raw))
        with self.assertRaisesRegex(ValueError, 'Unsupported manifest schema'):
            adapter.read(run)

    def test_modern_validation_and_discovery_retained(self):
        run = create_run(self.root / 'modern')
        self.assertEqual(adapter.discover(self.root), [str(run)])
        self.assertEqual(adapter.read(run)['schema_version'], 2)
        (run / 'tables/org2.tsv').write_text('modified')
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            adapter.read(run)
