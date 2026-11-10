import csv
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

from test_collect_run_summary import create_run, script

ROOT = Path(__file__).resolve().parents[1]
COMMAND = ['python3', str(ROOT / 'src/report/link_manifest_runs.py')]


def train_text(value):
    return (f'# substitution percent identity: {value}\n'
            '# substitution percent identity: 99\n'
            + ''.join('#last -' + option + ' 1\n' for option in ('t', 'a', 'A', 'b', 'B', 'S'))
            + '# score matrix (query letters = columns, reference letters = rows):\n'
            + ''.join(letter + ' 1 0 0 0\n' for letter in 'ACGT'))


def legacy(root, name='Out1_In2_In3', date='20251027', intermediate=False,
           identities=(90, 91, 97), genera=('Pan', 'Homo', 'Homo')):
    trio = root / name
    run = trio / date
    (run / 'metadata').mkdir(parents=True)
    names = name.split('_')
    manifest = {'organisms': [dict(slot=f'org{i}', short_name=short,
                                  raw_organism_name=f'{genus} species',
                                  metadata_json='/stale/absolute/path')
                              for i, (short, genus) in enumerate(zip(names, genera), 1)]}
    (run / 'metadata/metadata_manifest.json').write_text(json.dumps(manifest))
    target = run / 'intermediateFiles' if intermediate else run
    target.mkdir(exist_ok=True)
    for (ref, query), value in zip(((0, 1), (0, 2), (1, 2)), identities):
        (target / f'{names[ref]}2{names[query]}_{date}.train').write_text(train_text(value))
    return trio, run


def snapshot(root):
    if not root.exists():
        return None
    return {p.name: ('link', str(p.readlink())) if p.is_symlink()
            else ('file', p.read_bytes()) if p.is_file() else ('dir', None)
            for p in root.iterdir()}


class LinkManifestTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / 'tmp/codex', prefix='link-identity-')
        self.root = Path(self.temp.name)
        self.inputs = self.root / 'inputs'
        self.inputs.mkdir()
        self.output = self.root / 'links'

    def tearDown(self):
        self.temp.cleanup()

    def cli(self, *options, inputs=None, success=True):
        result = subprocess.run(COMMAND + ['--output-dir', str(self.output)]
                                + list(options) + [str(p) for p in (inputs or [self.inputs])],
                                capture_output=True, text=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def records(self):
        with (self.output / 'link_trio_filter.csv').open() as handle:
            return {row['trio_id']: row for row in csv.DictReader(handle)}

    def test_legacy_layout_latest_and_direct_input(self):
        trio, _ = legacy(self.inputs)
        _, latest = legacy(self.inputs, date='20261031101732', intermediate=True)
        (trio / 'bench_20270101').mkdir()
        (self.inputs / 'archive').mkdir()
        legacy(self.inputs / 'archive', name='Hidden1_Hidden2_Hidden3')
        self.cli(inputs=[trio, trio])
        self.assertEqual((self.output / trio.name).resolve(), trio)
        self.assertFalse((self.output / trio.name).readlink().is_absolute())
        self.assertEqual(self.records()[trio.name]['run_path'], str(latest))
        self.cli(inputs=[self.inputs])
        self.assertEqual(len(self.records()), 1)
        settings = json.loads((self.output / 'link_run_settings.json').read_text())
        self.assertEqual(settings['reused'], 1)
        self.assertEqual(settings['identity_source'], 'last-train-penultimate')
        self.assertEqual(self.records()[trio.name]['identity_source'], 'last-train-penultimate')

    def test_mixed_modern_legacy_and_latest_complete(self):
        trio, _ = legacy(self.inputs)
        old = create_run(self.inputs / 'modern', date='20250101')
        latest = create_run(self.inputs / 'modern', date='20260101')
        create_run(self.inputs / 'modern', date='20270101')
        manifest = script.runs.load(self.inputs / 'modern/20270101')
        manifest['status'] = 'pending'
        script.runs.atomic_json(self.inputs / 'modern/20270101' / script.runs.MANIFEST, manifest)
        self.cli()
        self.assertEqual((self.output / 'fixture-trio').resolve(), latest)
        self.assertEqual((self.output / trio.name).resolve(), trio)
        self.assertEqual(self.records()['fixture-trio']['layout'], 'v2')
        self.assertNotEqual(old, latest)

    def test_gap_waiver_and_absolute_bounds(self):
        retained, _ = legacy(self.inputs)
        excluded, _ = legacy(self.inputs, name='Alt1_Alt2_Alt3', genera=('Same',) * 3)
        self.cli('--min-idt-gap-between-out-in', '20', '--outgroup-idt-threshold', '89')
        rows = self.records()
        self.assertEqual(rows[retained.name]['min_idt_gap_between_out_in'], 'waived')
        self.assertEqual(rows[excluded.name]['included'], 'False')
        other = self.root / 'all-excluded'
        self.output = other
        self.cli('--min-idt-gap-between-out-in', '20', '--outgroup-idt-threshold', '90')
        self.assertFalse(any(p.is_symlink() for p in other.iterdir()))
        self.assertEqual(self.records()[retained.name]['outgroup_idt_threshold'], 'failed')

    def test_decimal_gap_and_common_override_order(self):
        trio, _ = legacy(self.inputs, identities=(91.8, 90, 95), genera=('Same',) * 3)
        self.cli('--min-idt-gap-between-out-in', '3.2')
        self.assertEqual(self.records()[trio.name]['included'], 'False')
        for i, args in enumerate([
            ['--idt-threshold', '99', '--ingroup-idt-threshold', '94', '--outgroup-idt-threshold', '89'],
            ['--outgroup-idt-threshold', '89', '--ingroup-idt-threshold', '94', '--idt-threshold', '99']]):
            self.output = self.root / f'ordered{i}'
            self.cli(*args, '--min-idt-gap-between-out-in', '3.19')
            self.assertTrue((self.output / trio.name).is_symlink())
            effective = json.loads((self.output / 'link_run_settings.json').read_text())['effective_thresholds']
            self.assertEqual(effective['ingroup_idt_threshold'], 94)

    def test_missing_identity_recorded_and_no_filter_disabled(self):
        trio, run = legacy(self.inputs)
        next(run.glob('*In22In3*.train')).unlink()
        self.cli('--ingroup-idt-threshold', '80')
        self.assertEqual(self.records()[trio.name]['ingroup_idt_threshold'], 'missing_identity')
        self.output = self.root / 'unfiltered'
        self.cli()
        self.assertTrue((self.output / trio.name).is_symlink())
        self.assertIn('missing training file', self.records()[trio.name]['reason'])

    def test_tighter_filter_preserves_all_existing_collection_links(self):
        first, _ = legacy(self.inputs)
        second, _ = legacy(self.inputs, name='Other1_Other2_Other3')
        self.cli()
        before = {p.name: p.readlink() for p in self.output.iterdir() if p.is_symlink()}
        result = self.cli('--ingroup-idt-threshold', '100')
        after = {p.name: p.readlink() for p in self.output.iterdir() if p.is_symlink()}
        self.assertEqual(before, after)
        self.assertEqual(set(after), {first.name, second.name})
        for row in self.records().values():
            self.assertEqual(row['included'], 'False')
            self.assertEqual(row['link_action'], 'kept_existing')
        settings = json.loads((self.output / 'link_run_settings.json').read_text())
        self.assertEqual(settings['created'], 0)
        self.assertEqual(settings['reused'], 0)
        self.assertEqual(settings['excluded'], 2)
        self.assertEqual(settings['preserved_existing'], 2)
        self.assertIn('preserved_existing: 2', result.stdout)
        self.cli('--ingroup-idt-threshold', '100')
        self.assertEqual(before, {
            p.name: p.readlink() for p in self.output.iterdir() if p.is_symlink()
        })

    def test_external_links_preserved_while_filtered_additions_are_created(self):
        trio, _ = legacy(self.inputs)
        excluded, _ = legacy(self.inputs, name='Low1_Low2_Low3', identities=(60, 65, 70))
        self.output.mkdir()
        unknown = self.output / 'unknown'
        unknown.symlink_to(trio)
        dangling = self.output / 'dangling'
        dangling.symlink_to('../missing-source')
        marker = self.output / 'keep.txt'
        marker.write_text('existing content')
        before = {p.name: p.readlink() for p in (unknown, dangling)}
        for expected_action in ('created', 'reused'):
            result = self.cli('--ingroup-idt-threshold', '80')
            self.assertTrue((self.output / trio.name).is_symlink())
            self.assertFalse((self.output / excluded.name).exists())
            self.assertEqual(before, {p.name: p.readlink() for p in (unknown, dangling)})
            self.assertEqual(marker.read_text(), 'existing content')
            rows = self.records()
            self.assertEqual(set(rows), {trio.name, excluded.name})
            self.assertEqual(rows[trio.name]['link_action'], expected_action)
            self.assertEqual(rows[excluded.name]['included'], 'False')
            self.assertEqual(rows[excluded.name]['link_action'], 'skipped')
            settings = json.loads((self.output / 'link_run_settings.json').read_text())
            self.assertEqual(settings['preserved_existing'], 2)
            self.assertEqual(settings['created'], int(expected_action == 'created'))
            self.assertEqual(settings['reused'], int(expected_action == 'reused'))
            self.assertIn('preserved_existing: 2', result.stdout)

    def test_wrong_target_and_directory_collisions_still_fail_before_writes(self):
        trio, _ = legacy(self.inputs)
        self.output.mkdir()
        target = self.output / trio.name
        target.symlink_to('../missing-source')
        before = snapshot(self.output)
        self.cli('--ingroup-idt-threshold', '100', success=False)
        self.assertEqual(before, snapshot(self.output))
        target.unlink()
        target.mkdir()
        before = snapshot(self.output)
        self.cli('--ingroup-idt-threshold', '100', success=False)
        self.assertEqual(before, snapshot(self.output))

    def test_collisions_duplicate_keys_and_equal_input_no_writes(self):
        trio, _ = legacy(self.inputs)
        self.output.mkdir()
        (self.output / trio.name).write_text('keep')
        before = snapshot(self.output)
        self.cli(success=False)
        self.assertEqual(before, snapshot(self.output))
        other = self.root / 'other'
        other.mkdir()
        legacy(other)
        self.cli(inputs=[self.inputs, other], success=False)
        self.assertEqual(before, snapshot(self.output))
        self.output = self.inputs
        before = snapshot(self.output)
        self.cli(success=False)
        self.assertEqual(before, snapshot(self.output))

    def test_output_subtree_and_cycle_ignored(self):
        trio, _ = legacy(self.inputs)
        self.output = self.inputs / 'collection'
        self.cli()
        (self.inputs / 'cycle').symlink_to(self.inputs)
        self.cli()
        self.assertEqual(set(self.records()), {trio.name})

    def test_invalid_manifest_and_ambiguous_train_no_writes(self):
        trio, run = legacy(self.inputs)
        (run / 'intermediateFiles').mkdir()
        train = next(run.glob('*.train'))
        (run / 'intermediateFiles' / train.name).write_bytes(train.read_bytes())
        self.cli(success=False)
        self.assertFalse(self.output.exists())
        (run / 'intermediateFiles' / train.name).unlink()
        path = run / 'metadata/metadata_manifest.json'
        original = json.loads(path.read_text())
        path.write_text(json.dumps(dict(original, schema_version=1)))
        self.cli(success=False)
        self.assertFalse(self.output.exists())
        path.write_text('not json')
        self.cli(success=False)
        self.assertFalse(self.output.exists())

    def test_invalid_options_and_missing_run(self):
        (self.inputs / 'Empty1_Empty2_Empty3').mkdir()
        self.cli(success=False)
        for value in ('-1', '101', 'nan', 'inf'):
            self.cli('--idt-threshold', value, success=False)
        self.assertFalse(self.output.exists())

    def test_modern_trio_named_parent_and_nonnumeric_run(self):
        run = create_run(self.inputs / 'Modern1_Modern2_Modern3', date='arbitrary-run')
        self.cli()
        self.assertEqual((self.output / 'fixture-trio').resolve(), run)
        self.assertEqual(self.records()['fixture-trio']['layout'], 'v2')

    def test_modern_pending_only_does_not_break_other_sources(self):
        run = create_run(self.inputs / 'Modern1_Modern2_Modern3', date='arbitrary-run')
        manifest = script.runs.load(run)
        manifest['status'] = 'pending'
        script.runs.atomic_json(run / script.runs.MANIFEST, manifest)
        trio, _ = legacy(self.inputs)
        self.cli()
        self.assertEqual(set(self.records()), {trio.name})

    def test_direct_modern_run_does_not_include_nested_runs(self):
        run = create_run(self.inputs / 'modern')
        nested = create_run(run / 'nested')
        manifest = script.runs.load(nested)
        manifest['trio_key'] = 'nested-trio'
        script.runs.atomic_json(nested / script.runs.MANIFEST, manifest)
        self.cli(inputs=[run])
        self.assertEqual(set(self.records()), {'fixture-trio'})

    def test_penultimate_identity_reverses_final_value_threshold_decision(self):
        trio, _ = legacy(self.inputs)  # analytical97, final99
        create_run(self.inputs / 'modern')  # same distinct analytical/final identities
        self.cli('--ingroup-idt-threshold', '98')
        rows = self.records()
        self.assertEqual(set(rows), {trio.name, 'fixture-trio'})
        for row in rows.values():
            self.assertEqual(row['idt_23'], '97.0')
            self.assertEqual(row['included'], 'False')
            self.assertEqual(row['ingroup_idt_threshold'], 'failed')
            self.assertEqual(row['identity_source'], 'last-train-penultimate')
