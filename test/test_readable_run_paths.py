import copy
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src/lib'))
import run_manifest as manifest


class ReadableRunPathsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / 'tmp/codex')
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.value = dict(schema_version=2, run_key='a' * 64, trio_key='b' * 64,
                          date='20260101', options={}, created_at='now',
                          artifacts=[], status='pending', organisms=[
            dict(slot=f'org{i}', accession=f'GCA_00000000{i}.{i}',
                 raw_organism_name=name, short_name='same',
                 fasta_path=f'/inputs/local_{i}/genomic.fa',
                 fasta_sha256=str(i) * 64, gff_sha256=None)
            for i, name in enumerate(('Alpha beta', 'Alpha berta', 'Gamma delta'), 1)])

    def path(self, value=None):
        return manifest.resolve_run(self.base, value or self.value)

    def save(self, run, value=None):
        manifest.prepare(run, copy.deepcopy(value or self.value))

    def test_full_names_versions_and_slot_order(self):
        run = self.path()
        self.assertEqual(run.parent.name,
            'Alpha_beta_GCA_000000001.1__Alpha_berta_GCA_000000002.2__Gamma_delta_GCA_000000003.3')
        self.assertEqual(run.name, '20260101-' + 'a' * 12)
        reversed_value = copy.deepcopy(self.value)
        reversed_value['organisms'].reverse()
        self.assertEqual(self.path(reversed_value), run)
        reversed_value['organisms'][0]['accession'] = 'GCA_000000003.4'
        self.assertNotEqual(self.path(reversed_value).parent, run.parent)

    def test_local_fallback_sanitization_and_equal_name_budget(self):
        for entry in self.value['organisms']:
            entry.pop('raw_organism_name')
            entry['accession'] = None
        self.value['organisms'][0]['ncbi_full_name'] = ' /é Alpha:::beta__ '
        self.value['organisms'][1]['fasta_path'] = '/inputs/!!!/x.fa'
        run = self.path()
        self.assertIn('Alpha_beta_local-' + '1' * 12, run.parent.name)
        self.assertIn('__org2_local-' + '2' * 12, run.parent.name)
        self.assertIn('__local_3_local-' + '3' * 12, run.parent.name)
        for entry in self.value['organisms']:
            entry['raw_organism_name'] = 'X' * 400
        with patch.object(manifest.os, 'pathconf', return_value=120):
            component = self.path().parent.name
        self.assertLessEqual(len(component.encode()), 120)
        self.assertEqual(len(set(len(p.split('_local-')[0]) for p in component.split('__'))), 1)

    def test_length_errors_preserve_tokens_and_label(self):
        with patch.object(manifest.os, 'pathconf', return_value=30):
            with self.assertRaisesRegex(ValueError, 'accession|tokens'):
                self.path()
        self.value['date'] = 'x' * 244
        with self.assertRaisesRegex(ValueError, 'label|length'):
            self.path()

    def test_old_and_new_incomplete_resume_ignore_display_and_input_paths(self):
        for old in (True, False):
            with self.subTest(old=old):
                run = (self.base / self.value['trio_key'] /
                       (self.value['date'] + '-' + self.value['run_key'])) if old else self.path()
                self.save(run)
                changed = copy.deepcopy(self.value)
                for entry in changed['organisms']:
                    entry['raw_organism_name'] = 'Renamed species'
                    entry['fasta_path'] = '/renamed/file.fa'
                self.assertEqual(self.path(changed), run)
                alias = self.base / 'alias'
                alias.symlink_to(run.parent, target_is_directory=True)
                self.assertEqual(self.path(changed), run)
                alias.unlink()
                shutil.rmtree(run.parent)

    def test_duplicate_matching_physical_runs_are_ambiguous(self):
        first = self.base / 'old' / 'run'
        second = self.base / 'new' / 'run'
        self.save(first)
        self.save(second)
        with self.assertRaisesRegex(ValueError, 'Ambiguous.*' + str(first)):
            self.path()

    def test_external_and_deeper_aliases_do_not_hijack_resume(self):
        output = self.base / 'output'
        output.mkdir()
        for target in (self.base / 'external' / 'run', output / 'nested' / 'deeper' / 'run'):
            self.save(target)
            alias = output / 'alias'
            alias.symlink_to(target.parent, target_is_directory=True)
            expected = manifest.resolve_run(output, self.value)
            self.assertEqual(expected.parent.parent, output)
            self.assertNotEqual(expected, target)
            alias.unlink()

    def test_out_of_scope_readable_parent_rejected_before_child_creation(self):
        output = self.base / 'output'
        output.mkdir()
        expected = manifest.resolve_run(output, self.value)
        external = self.base / 'external'
        external.mkdir()
        expected.parent.symlink_to(external, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'destination'):
            manifest.resolve_run(output, self.value)
        self.assertEqual(list(external.iterdir()), [])

    def test_same_depth_parent_alias_rejected_before_child_creation(self):
        expected = self.path()
        other = self.base / 'different-parent'
        other.mkdir()
        expected.parent.symlink_to(other, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, 'destination'):
            self.path()
        self.assertEqual(list(other.iterdir()), [])

    def test_short_collision_extends_and_exhaustion_fails(self):
        candidate = self.path()
        other = copy.deepcopy(self.value)
        other['run_key'] = 'a' * 12 + 'b' * 52
        self.save(candidate, other)
        self.assertEqual(self.path().name, '20260101-' + 'a' * 16)
        for size in range(16, 65, 4):
            self.save(candidate.with_name('20260101-' + 'a' * size), other)
        with self.assertRaisesRegex(ValueError, 'suffix|collision'):
            self.path()

    def test_unknown_and_corrupt_candidates_are_preserved(self):
        run = self.path()
        run.mkdir(parents=True)
        marker = run / 'unknown'
        marker.write_text('keep')
        with self.assertRaisesRegex(ValueError, 'manifest'):
            self.path()
        marker.unlink()
        (run / '.run.lock').touch()
        self.assertEqual(self.path(), run)
        (run / 'metadata').mkdir()
        (run / manifest.MANIFEST).write_text('{broken')
        with self.assertRaises(ValueError):
            self.path()
        self.assertEqual((run / manifest.MANIFEST).read_text(), '{broken')

    def test_unrelated_malformed_manifest_structures_do_not_block_allocation(self):
        for payload in ([], dict(self.value, organisms=[None])):
            path = self.base / 'unrelated' / 'broken' / manifest.MANIFEST
            manifest.atomic_json(path, payload)
            run = self.path()
            self.assertEqual(run.name, '20260101-' + 'a' * 12)
            manifest.atomic_json(run / manifest.MANIFEST, payload)
            with self.assertRaises(ValueError):
                self.path()
            self.assertEqual(json.loads((run / manifest.MANIFEST).read_text()), payload)
            shutil.rmtree(run.parent)

    def test_wrong_label_or_key_cannot_trigger_recovery(self):
        run = self.path()
        self.save(run)
        partial = run / 'intermediateFiles' / 'partial'
        partial.parent.mkdir()
        partial.write_text('keep')
        for field in ('date', 'run_key'):
            changed = copy.deepcopy(self.value)
            changed[field] = 'different'
            with self.assertRaises(ValueError):
                manifest.prepare(run, changed)
            self.assertEqual(partial.read_text(), 'keep')

    def test_explicit_destination_is_checked_before_initialization(self):
        wrong = self.base / 'arbitrary' / 'run'
        with self.assertRaisesRegex(ValueError, 'destination'):
            manifest.validate_destination(self.base, wrong, self.value)
        self.assertFalse(wrong.exists())
        manifest.validate_destination(self.base, self.path(), self.value)

    def test_changed_full_keys_and_labels_get_distinct_runs(self):
        first = self.path()
        self.save(first)
        changed = copy.deepcopy(self.value)
        changed['run_key'] = 'c' * 64
        self.assertEqual(self.path(changed).parent, first.parent)
        self.assertNotEqual(self.path(changed), first)
        changed = copy.deepcopy(self.value)
        changed['date'] = '20260102'
        self.assertNotEqual(self.path(changed), first)

    def test_mixed_discovery_keeps_latest_per_full_trio_key(self):
        for parent, key, stamp in (('oldhash', 'b' * 64, '2020'),
                                   ('readable', 'b' * 64, '2022'),
                                   ('sameaccessions', 'c' * 64, '2021')):
            value = copy.deepcopy(self.value)
            value.update(trio_key=key, created_at=stamp, status='complete')
            self.save(self.base / parent / 'run', value)
        self.assertEqual(set(manifest.discover(self.base)),
                         {self.base / 'readable/run', self.base / 'sameaccessions/run'})
