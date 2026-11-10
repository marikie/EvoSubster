"""Exercise real pipeline orchestration with bounded external-tool fixtures."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class StagedPipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='last-stages-')
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        binary = self.folder / 'bin'
        binary.mkdir()
        program = '''#!/usr/bin/env python3
import os, pathlib, sys, time
name = pathlib.Path(sys.argv[0]).name
if sys.argv[1:] == ['--version']:
    print('lastal ' + os.environ.get('FIXTURE_LAST_VERSION', '1654'))
    sys.exit(0)
folder = pathlib.Path(os.environ['FIXTURE_DIR'])
with (folder / 'events').open('a') as log: log.write(name + '\\n')
if name == 'lastdb':
    pathlib.Path(sys.argv[-2] + '.prj').write_text('db')
    sys.exit(0)
if name == 'last-train':
    print('# train')
    sys.exit(0)
if name == 'lastal':
    (folder / 'lastal-running').touch()
    print('##maf version=1\\n\\na score=4\\ns ref 0 4 + 4 ACGT\\ns query 0 4 + 4 ACGT\\n', flush=True)
    time.sleep(.1)
    (folder / 'lastal-running').unlink()
elif name == 'last-split':
    if (folder / 'lastal-running').exists():
        print('split stages overlap', file=sys.stderr)
        sys.exit(23)
    print(sys.stdin.read(), end='')
elif name == 'maf-linked':
    print(pathlib.Path(sys.argv[-1]).read_text(), end='')
if os.environ.get('SKIP_TOOL') == name:
    print('last-split: skipping sequence chr1 (memory)', file=sys.stderr)
if os.environ.get('FAIL_TOOL') == name:
    print('fixture failure diagnostic', file=sys.stderr)
    sys.exit(17)
'''
        for name in ('lastdb', 'last-train', 'lastal', 'last-split', 'maf-linked'):
            path = binary / name
            path.write_text(program)
            path.chmod(0o755)
        self.env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ['PATH'],
                        FIXTURE_DIR=str(self.folder), LASTAL_PIPELINE_MODE='staged', THREAD_NUM='4')
        for key in ('EVOSUBSTER_LAST_OPTIONS', 'LASTAL_THREADS', 'LASTDB_ALIGNMENT_SEED'):
            self.env.pop(key, None)
        for name in ('ref.fa', 'query.fa'):
            (self.folder / name).write_text('>seq\nACGT\n')
        (self.folder / 'train').write_text('# parameters\n')
        self.output = self.folder / 'output.maf'
        self.command = ['bash', str(ROOT / 'src/align/one2one.sh'), 'date',
                        str(self.folder / 'ref.fa'), str(self.folder / 'query.fa'),
                        str(self.folder / 'db'), str(self.folder / 'train'), str(self.output)]

    def run_pipeline(self, **env):
        return subprocess.run(self.command, env=dict(self.env, **env), capture_output=True, text=True)

    def events(self):
        return (self.folder / 'events').read_text().splitlines()

    def test_stages_do_not_overlap_and_record_resources(self):
        result = self.run_pipeline()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(self.output.is_file())
        self.assertTrue(Path(str(self.output) + '.complete.json').is_file())
        for stage in ('lastal', 'last-split', 'maf-linked'):
            timing = next(self.folder.glob('output.maf.' + stage + '.time.*'))
            self.assertIn('max_rss_kb=', timing.read_text())
            self.assertIn('exit_status=0', timing.read_text())

    def test_downstream_failure_preserves_completed_stage_and_resumes_without_db(self):
        result = self.run_pipeline(FAIL_TOOL='last-split')
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output.exists())
        embedded = Path(str(self.output) + '.lastal.maf')
        self.assertTrue(Path(str(embedded) + '.complete.json').exists())
        self.assertTrue(any('fixture failure diagnostic' in p.read_text()
                            for p in self.folder.glob('output.maf.last-split.stderr.*')))
        shutil.rmtree(self.folder / 'db.seed-RY128')
        result = self.run_pipeline()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.events().count('lastal'), 1)
        self.assertEqual(self.events().count('lastdb'), 1)

    def test_failed_or_skipping_stage_never_publishes_output(self):
        for stage in ('lastal', 'last-split', 'maf-linked'):
            with self.subTest(stage=stage):
                result = self.run_pipeline(FAIL_TOOL=stage)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.output.exists())
        for stage in ('lastal', 'last-split'):
            with self.subTest(skipping=stage):
                for path in self.folder.glob('output.maf*.complete.json'):
                    path.unlink()
                result = self.run_pipeline(SKIP_TOOL=stage)
                self.assertNotEqual(result.returncode, 0)
                self.assertFalse(self.output.exists())
                self.assertIn('incomplete alignments are not accepted', result.stderr)

    def test_changed_sealed_stage_is_rerun(self):
        self.run_pipeline(FAIL_TOOL='last-split')
        Path(str(self.output) + '.lastal.maf').write_text('partial')
        result = self.run_pipeline()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.events().count('lastal'), 2)
        self.assertNotIn('partial', self.output.read_text())

    def test_unsealed_final_is_not_left_published_when_retry_fails(self):
        self.output.write_text('partial final output')
        result = self.run_pipeline(FAIL_TOOL='lastal')
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output.exists())

    def test_changed_input_invalidates_checkpoint(self):
        self.run_pipeline(FAIL_TOOL='last-split')
        (self.folder / 'query.fa').write_text('>query\nACTT\n')
        result = self.run_pipeline()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.events().count('lastal'), 2)

    def test_last_version_change_invalidates_checkpoint_and_final(self):
        result = self.run_pipeline(FAIL_TOOL='last-split', FIXTURE_LAST_VERSION='1648')
        self.assertNotEqual(result.returncode, 0)
        result = self.run_pipeline(FIXTURE_LAST_VERSION='1654')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.events().count('lastal'), 2)
        record = json.loads(Path(str(self.output) + '.complete.json').read_text())
        self.assertEqual(record['lastal_version'], 'lastal 1654')
        result = self.run_pipeline(FIXTURE_LAST_VERSION='1655')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.events().count('lastal'), 3)

    def test_moved_checkpoint_reuses_logical_inputs_and_new_thread_count(self):
        result = self.run_pipeline(FAIL_TOOL='last-split')
        self.assertNotEqual(result.returncode, 0)
        moved = self.folder / 'moved'
        moved.mkdir()
        for name in ('ref.fa', 'query.fa', 'train', 'output.maf.lastal.maf',
                     'output.maf.lastal.maf.complete.json', 'output.maf.last_options.json'):
            shutil.move(self.folder / name, moved / name)
        for index in (3, 4, 6, 7):
            self.command[index] = str(moved / Path(self.command[index]).name)
        result = self.run_pipeline(LASTAL_THREADS='1')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.events().count('lastal'), 1)
        seal = json.loads((moved / 'output.maf.complete.json').read_text())
        self.assertEqual([item['role'] for item in seal['inputs']], ['reference', 'query', 'train'])
        self.assertNotIn(str(self.folder), json.dumps(seal))

    def test_legacy_seal_still_requires_exact_settings_and_input_paths(self):
        result = self.run_pipeline(FAIL_TOOL='last-split')
        self.assertNotEqual(result.returncode, 0)
        checkpoint = Path(str(self.output) + '.lastal.maf')
        helper = ROOT / 'src/lib/last_options.py'
        import importlib.util
        spec = importlib.util.spec_from_file_location('legacy_last_options', helper)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        from unittest.mock import patch
        with patch.dict(os.environ, self.env, clear=True):
            seal = dict(schema=1, settings=module.config(), lastal_version='lastal 1654',
                        sha256=module.file_digest(checkpoint),
                        inputs={str((self.folder / name).resolve()): module.file_digest(self.folder / name)
                                for name in ('ref.fa', 'query.fa', 'train')})
        Path(str(checkpoint) + '.complete.json').write_text(json.dumps(seal))
        verify = ['python3', str(helper), 'verify-output', str(checkpoint)]
        for name in ('ref.fa', 'query.fa', 'train'):
            verify.extend(['--input', str(self.folder / name)])
        changed = subprocess.run(verify, env=dict(self.env, LASTAL_THREADS='1'),
                                 capture_output=True, text=True)
        self.assertNotEqual(changed.returncode, 0)
        moved_ref = self.folder / 'moved-ref.fa'
        shutil.copyfile(self.folder / 'ref.fa', moved_ref)
        moved = list(verify)
        moved[5] = str(moved_ref)
        changed = subprocess.run(moved, env=self.env, capture_output=True, text=True)
        self.assertNotEqual(changed.returncode, 0)
        result = self.run_pipeline()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.events().count('lastal'), 1)

    def test_streaming_publishes_only_successful_final_seal(self):
        # The fixture's overlap rejection is only relevant to staged execution.
        split = self.folder / 'bin/last-split'
        split.write_text(split.read_text().replace("if (folder / 'lastal-running').exists():", "if False:"))
        result = self.run_pipeline(LASTAL_PIPELINE_MODE='streaming')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(Path(str(self.output) + '.complete.json').is_file())
        self.assertFalse(Path(str(self.output) + '.raw.complete.json').exists())

    def test_final_resume_checks_seal_and_settings_before_building_database(self):
        result = self.run_pipeline()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        shutil.rmtree(self.folder / 'db.seed-RY128')
        result = self.run_pipeline()
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.events().count('lastdb'), 1)
        result = self.run_pipeline(LASTAL_THREADS='1')
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.events().count('lastal'), 1)


if __name__ == '__main__':
    unittest.main()
