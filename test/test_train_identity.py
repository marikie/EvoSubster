import json
import subprocess
from pathlib import Path
import sys
import tempfile
import unittest

from test_collect_run_summary import create_run, script
from test_link_manifest_runs import train_text

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src/lib'))
from train_identity import IDENTITY_SOURCE, parse_identity
import run_manifest as runs


class TrainIdentityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / 'tmp/codex', prefix='penultimate-')
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_penultimate_distinct_from_final(self):
        path = self.root / 'pair.train'
        path.write_text('# substitution percent identity: 15\n' + train_text(84.5))
        self.assertEqual(parse_identity(path), 84.5)
        self.assertIn('# substitution percent identity: 99', path.read_text())

    def test_malformed_penultimate_never_falls_back(self):
        path = self.root / 'pair.train'
        for value in ('NaN', '101', '-1', 'bad', ''):
            path.write_text('# substitution percent identity: 85\n' + train_text(90).replace(
                '# substitution percent identity: 90', '# substitution percent identity: ' + value))
            with self.assertRaisesRegex(ValueError, 'Invalid penultimate'):
                parse_identity(path)

    def test_one_marker_and_incomplete_final_block_rejected(self):
        path = self.root / 'pair.train'
        path.write_text(train_text(90).replace('# substitution percent identity: 90\n', ''))
        with self.assertRaisesRegex(ValueError, 'At least two'):
            parse_identity(path)
        for text in (train_text(90).replace('#last -S 1\n',''),
                     train_text(90).replace('T 1 0 0 0\n','')):
            path.write_text(text)
            with self.assertRaisesRegex(ValueError, 'Incomplete LAST'):
                parse_identity(path)

    def old_run(self):
        run = create_run(self.root / 'old')
        manifest = runs.load(run)
        manifest['options']['historical_choice'] = 'preserve'
        train = next(a for a in manifest['artifacts'] if a['kind'] == 'train')
        manifest['selection_provenance'] = dict(path=train['path'], sha256=train['sha256'])
        for artifact in manifest['artifacts']:
            if artifact['kind'] == 'train':
                artifact.pop('identity_source')
                artifact['identity_percent'] = 99
        runs.atomic_json(run / runs.MANIFEST, manifest)
        return run, manifest

    def test_explicit_migration_preserves_every_other_field_and_is_idempotent(self):
        run, original = self.old_run()
        path = run / runs.MANIFEST
        before = path.read_bytes()
        with self.assertRaisesRegex(ValueError, 'migrate-identity'):
            runs.validated(run, complete=True)
        self.assertEqual(path.read_bytes(), before)
        result = subprocess.run([
            sys.executable, str(ROOT / 'src/lib/run_manifest.py'), 'migrate-identity', str(run)
        ], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        migrated = runs.validated(run, complete=True)
        backups = list(path.parent.glob(path.name + '.before-penultimate-*.bak'))
        self.assertEqual(len(backups), 1)
        self.assertEqual(backups[0].read_bytes(), before)
        for old, new in zip(original['artifacts'], migrated['artifacts']):
            if old['kind'] == 'train':
                self.assertEqual(new['identity_source'], IDENTITY_SOURCE)
                self.assertNotEqual(old['identity_percent'], new['identity_percent'])
                new.pop('identity_source')
                new['identity_percent'] = old['identity_percent']
        self.assertEqual(migrated, original)
        after = path.read_bytes()
        runs.migrate_identity(run)
        self.assertEqual(path.read_bytes(), after)
        self.assertEqual(len(list(path.parent.glob('*.bak'))), 1)

    def test_migration_refuses_bad_hash_and_truncated_training_without_writes(self):
        run, manifest = self.old_run()
        path = run / runs.MANIFEST
        train = next(a for a in manifest['artifacts'] if a['kind'] == 'train')
        file = run / train['path']
        file.write_text(file.read_text().replace('T 1 0 0 0\n', ''))
        before = path.read_bytes()
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            runs.migrate_identity(run)
        self.assertEqual(path.read_bytes(), before)
        self.assertFalse(list(path.parent.glob('*.bak')))
        train['sha256'] = runs.sha256(file)
        manifest['selection_provenance']['sha256'] = train['sha256']
        runs.atomic_json(path, manifest)
        before = path.read_bytes()
        with self.assertRaisesRegex(ValueError, 'Incomplete LAST'):
            runs.migrate_identity(run)
        self.assertEqual(path.read_bytes(), before)
        self.assertFalse(list(path.parent.glob('*.bak')))

    def test_migration_refuses_locked_run(self):
        import fcntl
        run, _ = self.old_run()
        path = run / runs.MANIFEST
        before = path.read_bytes()
        with (run / '.run.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(ValueError, 'already active'):
                runs.migrate_identity(run)
        self.assertEqual(path.read_bytes(), before)
        self.assertFalse(list(path.parent.glob('*.bak')))

    def test_identity_only_run_migrates_without_status_or_history_changes(self):
        run, manifest = self.old_run()
        manifest['status'] = 'identity_complete'
        manifest['artifacts'] = [a for a in manifest['artifacts'] if a['kind'] == 'train']
        runs.atomic_json(run / runs.MANIFEST, manifest)
        runs.migrate_identity(run)
        value = runs.validated(run)
        self.assertEqual(value['status'], 'identity_complete')
        self.assertEqual(value['options'], manifest['options'])
        self.assertEqual(value['selection_provenance'], manifest['selection_provenance'])
        self.assertTrue(all(a['identity_source'] == IDENTITY_SOURCE for a in value['artifacts']))

    def test_migration_rejects_one_marker_and_invalid_penultimate(self):
        for mode in ('one_marker', 'invalid_penultimate'):
            run, manifest = self.old_run()
            path = run / runs.MANIFEST
            train = next(a for a in manifest['artifacts'] if a['kind'] == 'train')
            file = run / train['path']
            text = file.read_text()
            marker = text.splitlines()[0] + '\n'
            replacement = '' if mode == 'one_marker' else '# substitution percent identity: bad\n'
            file.write_text(text.replace(marker, replacement, 1))
            train['sha256'] = runs.sha256(file)
            manifest['selection_provenance']['sha256'] = train['sha256']
            runs.atomic_json(path, manifest)
            before = path.read_bytes()
            with self.assertRaisesRegex(ValueError, 'At least two|Invalid penultimate'):
                runs.migrate_identity(run)
            self.assertEqual(path.read_bytes(), before)
            self.assertFalse(list(path.parent.glob('*.bak')))

    def test_migration_refuses_incorrect_old_stored_identity(self):
        run, manifest = self.old_run()
        train = next(a for a in manifest['artifacts'] if a['kind'] == 'train')
        train['identity_percent'] = 42
        path = run / runs.MANIFEST
        runs.atomic_json(path, manifest)
        before = path.read_bytes()
        with self.assertRaisesRegex(ValueError, 'requires verification/migration'):
            runs.migrate_identity(run)
        self.assertEqual(path.read_bytes(), before)
        self.assertFalse(list(path.parent.glob('*.bak')))
