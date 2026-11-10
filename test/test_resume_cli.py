import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
ENV_NAMES = ('LASTAL_THREADS', 'LASTAL_PIPELINE_MODE', 'EVOSUBSTER_REUSE_FROM',
             'EVOSUBSTER_REUSE_LAST_VERSION', 'EVOSUBSTER_REUSE_PLAN')


class ResumeCliTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='resume-cli-')
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.source = self.folder / 'old results'
        self.source.mkdir()
        self.env = os.environ.copy()
        for name in (*ENV_NAMES, 'EVOSUBSTER_LAST_OPTIONS', 'THREAD_NUM'):
            self.env.pop(name, None)

    def parse(self, args, env=None, deferred=False):
        script = '''source "$1"
shift
last_options_parse_cli "$@" || exit $?
python3 - "${LAST_REMAINING_ARGS[@]}" <<'END'
import json, os, sys
names = ('LASTAL_THREADS', 'LASTAL_PIPELINE_MODE', 'EVOSUBSTER_REUSE_FROM',
         'EVOSUBSTER_REUSE_LAST_VERSION', 'EVOSUBSTER_REUSE_PLAN')
print(json.dumps({'env': {name: os.environ.get(name) for name in names},
                  'remaining': sys.argv[1:],
                  'options': json.loads(os.environ['EVOSUBSTER_LAST_OPTIONS'])}))
END
'''
        effective = dict(self.env, **(env or {}))
        if deferred:
            effective['LAST_OPTIONS_DEFER_VALIDATE'] = '1'
        return subprocess.run(['bash', '-c', script, 'resume-cli',
                               str(ROOT / 'src/lib/last_options.sh'), *args],
                              cwd=self.folder, env=effective, capture_output=True, text=True)

    def test_cli_overrides_environment_and_preserves_existing_arguments(self):
        result = self.parse([
            'run-label', '--lastal-threads', '2', '--alignment-mode=staged',
            '--reuse-from', 'old results', '--reuse-last-version=allow-existing-pairs',
            '--reuse-plan', '--lastal-options=-H2', '--out-dir', 'new results'],
            env={'LASTAL_THREADS': '8', 'LASTAL_PIPELINE_MODE': 'streaming',
                 'EVOSUBSTER_REUSE_FROM': '/does/not/exist',
                 'EVOSUBSTER_REUSE_LAST_VERSION': 'same', 'EVOSUBSTER_REUSE_PLAN': '0'})
        self.assertEqual(result.returncode, 0, result.stderr)
        parsed = json.loads(result.stdout)
        self.assertEqual(parsed['env'], {
            'LASTAL_THREADS': '2', 'LASTAL_PIPELINE_MODE': 'staged',
            'EVOSUBSTER_REUSE_FROM': str(self.source),
            'EVOSUBSTER_REUSE_LAST_VERSION': 'allow-existing-pairs',
            'EVOSUBSTER_REUSE_PLAN': '1'})
        self.assertEqual(parsed['remaining'], ['run-label', '--out-dir', 'new results'])
        self.assertEqual(parsed['options'], {'lastal': '-H2'})

    def test_all_value_options_accept_equals_and_separate_forms(self):
        for equals in (False, True):
            values = [('--lastal-threads', '3'), ('--alignment-mode', 'streaming'),
                      ('--reuse-from', 'old results'), ('--reuse-last-version', 'same')]
            args = [part for name, value in values
                    for part in ([name + '=' + value] if equals else [name, value])]
            with self.subTest(equals=equals):
                result = self.parse(args)
                self.assertEqual(result.returncode, 0, result.stderr)
                parsed = json.loads(result.stdout)['env']
                self.assertEqual(parsed['LASTAL_THREADS'], '3')
                self.assertEqual(parsed['EVOSUBSTER_REUSE_FROM'], str(self.source))

    def test_defaults_and_environment_fallback(self):
        result = self.parse([])
        self.assertEqual(result.returncode, 0, result.stderr)
        parsed = json.loads(result.stdout)['env']
        self.assertEqual(parsed['LASTAL_PIPELINE_MODE'], 'streaming')
        self.assertEqual(parsed['EVOSUBSTER_REUSE_LAST_VERSION'], 'same')
        self.assertIsNone(parsed['LASTAL_THREADS'])
        result = self.parse([], env={'LASTAL_THREADS': '4', 'LASTAL_PIPELINE_MODE': 'staged',
                                   'EVOSUBSTER_REUSE_FROM': 'old results',
                                   'EVOSUBSTER_REUSE_LAST_VERSION': 'allow-existing-pairs',
                                   'EVOSUBSTER_REUSE_PLAN': '1'})
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)['env']['EVOSUBSTER_REUSE_FROM'], str(self.source))

    def test_duplicates_rejected_across_argument_forms(self):
        for name, value in [('--lastal-threads', '2'), ('--alignment-mode', 'staged'),
                            ('--reuse-from', 'old results'), ('--reuse-last-version', 'same')]:
            with self.subTest(name=name):
                result = self.parse([name, value, name + '=' + value])
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(name, result.stderr)
        result = self.parse(['--reuse-plan', '--reuse-plan'])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('--reuse-plan', result.stderr)

    def test_missing_and_invalid_values_fail_even_when_validation_is_deferred(self):
        for args in [
            ['--lastal-threads'], ['--lastal-threads='], ['--lastal-threads', '0'],
            ['--lastal-threads', '-1'], ['--lastal-threads', '1.5'],
            ['--lastal-threads', 'abc'], ['--alignment-mode'], ['--alignment-mode='],
            ['--alignment-mode', 'unknown'], ['--reuse-last-version'],
            ['--reuse-last-version='], ['--reuse-last-version', 'any'],
            ['--reuse-from'], ['--reuse-from='], ['--reuse-from', '--reuse-plan'],
            ['--reuse-from', 'missing'], ['--reuse-plan=1'],
        ]:
            with self.subTest(args=args):
                result = self.parse(args, deferred=True)
                self.assertNotEqual(result.returncode, 0, result.stdout)
                self.assertIn(args[0].split('=')[0], result.stderr)
                self.assertFalse((self.folder / 'missing').exists())

    def test_invalid_environment_fallback_is_rejected(self):
        for name, value in [('LASTAL_THREADS', '0'), ('LASTAL_PIPELINE_MODE', 'other'),
                            ('EVOSUBSTER_REUSE_LAST_VERSION', 'any'),
                            ('EVOSUBSTER_REUSE_FROM', 'missing'), ('EVOSUBSTER_REUSE_PLAN', 'maybe')]:
            with self.subTest(name=name):
                result = self.parse([], env={name: value}, deferred=True)
                self.assertNotEqual(result.returncode, 0)

    def test_reuse_source_must_be_directory_and_is_never_created(self):
        plain_file = self.folder / 'plain-file'
        plain_file.write_text('not a result directory')
        result = self.parse(['--reuse-from', str(plain_file)])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('--reuse-from', result.stderr)
        self.assertEqual(sorted(path.name for path in self.folder.iterdir()), ['old results', 'plain-file'])

    def test_double_dash_preserves_positional_arguments(self):
        result = self.parse(['--', '--reuse-plan', '--alignment-mode=invalid'])
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)['remaining'],
                         ['--', '--reuse-plan', '--alignment-mode=invalid'])

    def test_both_entrypoints_reject_invalid_settings_before_output_creation(self):
        for entrypoint in ('sbst.sh', 'sbst_fromDwl.sh'):
            for args in [ ['--lastal-threads=0'], ['--alignment-mode=other'],
                          ['--reuse-last-version=any'], ['--reuse-from', 'missing'] ]:
                with self.subTest(entrypoint=entrypoint, args=args):
                    result = subprocess.run(['bash', str(ROOT / 'src' / entrypoint),
                                             '--out-dir', str(self.folder / 'new'), *args],
                                            cwd=self.folder, env=self.env,
                                            capture_output=True, text=True)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(args[0].split('=')[0], result.stderr)
                    self.assertFalse((self.folder / 'new').exists())
                    self.assertNotIn('Downloading', result.stdout)


if __name__ == '__main__':
    unittest.main()
