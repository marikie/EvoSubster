import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('last_options', ROOT / 'src/lib/last_options.py')
options = importlib.util.module_from_spec(spec)
spec.loader.exec_module(options)


class LastOptionsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='last-options-')
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.env_patch = patch.dict(os.environ, {}, clear=False)
        self.env_patch.start()
        self.addCleanup(self.env_patch.stop)
        for key in ('EVOSUBSTER_LAST_OPTIONS', 'THREAD_NUM', 'LASTDB_ALIGNMENT_SEED',
                    'LASTAL_QUERY_BATCH_SIZE', 'LASTAL_SPLIT_MEMORY'):
            os.environ.pop(key, None)

    def custom(self, **values):
        os.environ['EVOSUBSTER_LAST_OPTIONS'] = json.dumps(values)

    def test_quoted_values_and_overrides(self):
        self.custom(lastdb='-uRY4 -w 2', lastal='-H2 -C1 --split-b=32G',
                    **{'last-train': '--fixmat "matrix with spaces" --sample-number=1000',
                       'last-split': '--mismap=0.1'})
        config = options.config()
        self.assertEqual(config['lastdb_align'], ['-P', '8', '-c', '-u', 'RY4', '-w', '2'])
        self.assertEqual(config['lastal'].count('-H'), 1)
        self.assertEqual(config['lastal'][config['lastal'].index('-H') + 1], '2')
        self.assertNotIn('8T', config['lastal'])
        self.assertIn('matrix with spaces', config['last_train'])
        self.assertEqual(config['last_split'], ['-r', '--mismap', '0.1'])

    def test_rejects_operands_unknown_options_and_managed_settings(self):
        for tool, value in [
            ('lastdb', '-P4'), ('lastdb', '-vcP4'), ('lastdb', '-p'),
            ('lastal', '-p profiles.train'), ('lastal', '-fMAF'), ('lastal', '--split-f=TAB'),
            ('lastal', '--splice'), ('last-train', '-P 2'), ('last-split', '-r'),
            ('last-split', '--format=MAF+'), ('last-split', '-n'),
            ('lastal', 'input.fa'), ('lastal', '--'), ('lastal', '--unknown'),
            ('lastal', '-H'), ('lastal', '--split-b='), ('lastal', "-H '"),
            ('lastal', '--help'), ('last-train', '--codon'),
        ]:
            with self.subTest(tool=tool, value=value), self.assertRaises(ValueError):
                options.parse_options(tool, value)

    def test_environment_defaults_and_cli_precedence(self):
        os.environ.update(THREAD_NUM='3', LASTDB_ALIGNMENT_SEED='RY4',
                          LASTAL_QUERY_BATCH_SIZE='2M', LASTAL_SPLIT_MEMORY='16G')
        self.assertIn('2M', options.config()['lastal'])
        self.custom(lastdb='-uRY128', lastal='-i4M --split-b=32G')
        self.assertIn('RY128', options.config()['lastdb_align'])
        self.assertNotIn('2M', options.config()['lastal'])
        self.assertIn('4M', options.config()['lastal'])
        self.assertIn('32G', options.config()['lastal'])

    def test_seed_paths_cannot_become_database_directory_paths(self):
        self.custom(lastdb='-u "/tmp/seed definitions/RY4"')
        label = options.database_label()
        self.assertTrue(label.startswith('custom-'))
        self.assertNotIn('/', label)
        self.assertNotIn(' ', label)
        os.environ['LASTDB_ALIGNMENT_SEED'] = ''
        self.custom(lastdb='-uRY4')
        self.assertEqual(options.database_label(), 'RY4')

    def test_run_guard_and_training_cache_settings(self):
        output = self.folder / 'output'
        cache = self.folder / 'cache'
        self.custom(lastal='-H2')
        options.guard(output, 'run')
        options.guard(cache, 'training')
        options.guard(output, 'run')
        self.custom(lastal='-H3')
        options.guard(cache, 'training')  # Alignment settings do not alter training.
        with self.assertRaisesRegex(ValueError, 'settings differ'):
            options.guard(output, 'run')
        self.custom(**{'last-train': '--sample-number=1000'})
        with self.assertRaisesRegex(ValueError, 'settings differ'):
            options.guard(cache, 'training')

    def test_legacy_output_and_cache_not_reused_with_custom_options(self):
        for stage in ('run', 'training', 'database'):
            with self.subTest(stage=stage):
                folder = self.folder / stage
                folder.mkdir()
                original = folder / 'existing.train'
                original.write_text('keep existing data')
                self.custom(lastdb='-w2')
                with self.assertRaisesRegex(ValueError, 'no LAST settings record'):
                    options.guard(folder, stage)
                self.assertEqual(original.read_text(), 'keep existing data')
        legacy = self.folder / 'default'
        legacy.mkdir()
        (legacy / 'old.maf').write_text('legacy')
        self.custom()
        options.guard(legacy, 'run')

    def test_cli_validation_and_duplicate_detection_before_download(self):
        for script in ('sbst.sh', 'sbst_fromDwl.sh'):
            for args in [
                ['--lastal-options', "-H '"],
                ['--lastal-options=-H2', '--lastal-options=-H3'],
                ['--lastdb-options', '-P4'], ['--lastal-options'],
            ]:
                with self.subTest(script=script, args=args):
                    result = subprocess.run(['rtk', 'proxy', 'bash', str(ROOT / 'src' / script), *args],
                                            capture_output=True, text=True)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertNotIn('Downloading', result.stdout)
            result = subprocess.run(['rtk', 'proxy', 'bash', str(ROOT / 'src' / script),
                                     '--lastal-options=', '--help'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('--last-split-options', result.stdout)

    def make_tools(self):
        binary = self.folder / 'bin'
        binary.mkdir()
        program = '''#!/usr/bin/env python3
import json, os, pathlib, sys
tool = pathlib.Path(sys.argv[0]).name
with open(os.environ['OPTIONS_LOG'], 'a') as log:
    log.write(json.dumps([tool, sys.argv[1:]]) + '\\n')
if tool == 'lastdb':
    pathlib.Path(sys.argv[-2] + '.prj').touch()
elif tool == 'lastal':
    print('##maf version=1\\n\\na\\ns ref 0 3 + 3 ACG\\ns query 0 3 + 3 ACG\\n')
elif tool == 'last-train':
    print('# parameters')
elif tool == 'last-split':
    print(sys.stdin.read(), end='')
elif tool == 'maf-linked':
    print(pathlib.Path(sys.argv[-1]).read_text(), end='')
'''
        for name in ('lastdb', 'lastal', 'last-train', 'last-split', 'maf-linked'):
            path = binary / name
            path.write_text(program)
            path.chmod(0o755)
        os.environ['PATH'] = str(binary) + os.pathsep + os.environ['PATH']
        os.environ['OPTIONS_LOG'] = str(self.folder / 'arguments.jsonl')

    def test_actual_helper_argv_and_no_shell_execution(self):
        self.make_tools()
        payload = '$(touch {})'.format(self.folder / 'should-not-exist')
        self.custom(lastdb='-uRY4 -w2', lastal='-H2 -C1',
                    **{'last-train': '--fixmat ' + repr(payload), 'last-split': '-m0.1'})
        ref, query = self.folder / 'reference.fa', self.folder / 'query.fa'
        ref.write_text('>ref\nACG\n')
        query.write_text('>query\nACG\n')
        output = self.folder / 'fresh/output.maf'
        command = ['rtk', 'proxy', 'bash', str(ROOT / 'src/align/one2one.sh'), 'date',
                   str(ref), str(query), str(self.folder / 'db'),
                   str(self.folder / 'params.train'), str(output)]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        logged = {name: argv for name, argv in map(json.loads, Path(os.environ['OPTIONS_LOG']).read_text().splitlines())}
        self.assertIn('-uRY4', logged['lastdb'])
        self.assertIn('--fixmat=' + payload, logged['last-train'])
        self.assertIn('-H', logged['lastal'])
        self.assertEqual(logged['lastal'][logged['lastal'].index('-H') + 1], '2')
        self.assertEqual(logged['last-split'], ['-r', '-m', '0.1'])
        self.assertFalse((self.folder / 'should-not-exist').exists())
        self.assertTrue(output.is_file())
        self.custom(lastal='-H3')
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('settings differ', result.stderr)


if __name__ == '__main__':
    unittest.main()
