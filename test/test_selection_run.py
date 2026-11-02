import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

LIB = Path(__file__).resolve().parents[1] / 'src' / 'lib'
sys.path.insert(0, str(LIB))
import selection_run


class SelectionRunTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='trio settings ')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.tsv = self.root / 'selected_trios.tsv'
        self.tsv.write_text('in2_acc\tnote\tout_acc\tin1_acc\n'
                            'GCF_3.1\tany details\tGCA_1.1\tGCA_2.1\n')
        self.env = {k: v for k, v in os.environ.items()
                    if k not in selection_run.ENV_DEFAULTS and
                    k not in ('EVOSUBSTER_LAST_OPTIONS', 'THREAD_NUM')}

    def cli(self, *args, env=None):
        return subprocess.run([sys.executable, str(LIB / 'selection_run.py'), *args],
                              env=env or self.env, capture_output=True)

    def save(self):
        env = dict(self.env, LASTAL_QUERY_BATCH_SIZE='8K', LASTAL_SPLIT_MEMORY='94G',
                   EVOSUBSTER_LAST_OPTIONS=json.dumps({'lastal': '-H2', 'lastdb': '-uRY4'}))
        result = self.cli('save', str(self.root / 'selection_run.json'),
                          '--date', '20260818', '--threads', '1',
                          '--genome-dir', str(self.root / 'genomes'),
                          '--out-dir', str(self.root / 'out'),
                          '--train-cache-dir', str(self.root / 'cache'), env=env)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_rows_reordered_with_extra_columns_and_blank_lines(self):
        self.assertEqual(selection_run.read_trios(self.tsv), [['GCA_1.1', 'GCA_2.1', 'GCF_3.1']])
        self.tsv.write_text('out_acc\tin1_acc\tin2_acc\n\nGCA_1.1\tGCA_2.1\tGCF_3.1\n\n')
        self.assertEqual(len(selection_run.read_trios(self.tsv)), 1)

    def test_empty_inputs(self):
        for content in ('', '\n', 'out_acc\tin1_acc\tin2_acc\n'):
            self.tsv.write_text(content)
            self.assertEqual(selection_run.read_trios(self.tsv), [])

    def test_invalid_rows_fail_before_emitting_any_trio(self):
        for row in ('GCA_1.1\tGCA_1.1\tGCF_3.1',
                    'GCA_1\tGCA_2.1\tGCF_3.1', 'GCA_1.1\tGCA_2.1'):
            self.tsv.write_text('out_acc\tin1_acc\tin2_acc\n'
                                'GCA_1.1\tGCA_2.1\tGCF_3.1\n' + row + '\n')
            result = self.cli('trios', str(self.tsv))
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(result.stdout, b'')
            self.assertIn(b'line 3', result.stderr)

    def test_missing_and_duplicated_columns(self):
        for header in ('out_acc\tin1_acc', 'out_acc\tin1_acc\tin2_acc\tout_acc'):
            self.tsv.write_text(header + '\n')
            with self.assertRaisesRegex(ValueError, 'column'):
                selection_run.read_trios(self.tsv)

    def test_round_trip_restores_settings_over_ambient_environment(self):
        self.save()
        env = dict(self.env, LASTAL_QUERY_BATCH_SIZE='1M',
                   EVOSUBSTER_LAST_OPTIONS='{"lastal":"-H3"}')
        result = self.cli('resolve', '--trios', str(self.tsv), env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        values = result.stdout.decode().split('\0')
        self.assertEqual(values[:5], ['20260818', str(self.root / 'genomes'),
                                     str(self.root / 'out'), str(self.root / 'cache'), '1'])
        self.assertEqual(json.loads(values[5])['lastal'], '-H2')
        self.assertEqual(values[6:9], ['RY128', '8K', '94G'])

    def test_explicit_overrides_and_relocated_config(self):
        self.save()
        moved = self.root / 'moved.json'
        (self.root / 'selection_run.json').rename(moved)
        result = self.cli('resolve', '--trios', str(self.tsv), '--config', str(moved),
                          '--date', 'new', '--threads', '4', '--out-dir', '/tmp/new',
                          '--explicit-last', 'lastal',
                          env=dict(self.env, EVOSUBSTER_LAST_OPTIONS='{"lastal":""}'))
        self.assertEqual(result.returncode, 0, result.stderr)
        values = result.stdout.decode().split('\0')
        self.assertEqual((values[0], values[2], values[4]), ('new', '/tmp/new', '4'))
        self.assertEqual(json.loads(values[5])['lastal'], '')
        self.assertEqual(json.loads(values[5])['lastdb'], '-uRY4')

    def test_legacy_requires_date(self):
        result = self.cli('resolve', '--trios', str(self.tsv))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b'DATE is required', result.stderr)
        result = self.cli('resolve', '--trios', str(self.tsv), '--date', 'old')
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_malformed_config_never_falls_back(self):
        for value in ('{broken', '{"schema":99}', '{"schema":1}'):
            (self.root / 'selection_run.json').write_text(value)
            result = self.cli('resolve', '--trios', str(self.tsv), '--date', 'old')
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(result.stdout, b'')


if __name__ == '__main__':
    unittest.main()
