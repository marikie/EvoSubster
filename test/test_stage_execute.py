"""A missing or non-GNU time command cannot prevent a stage from running."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class StageTimingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='last-stage-timing-')
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)

    def run_stage(self, timer, status=0):
        command = [sys.executable, str(ROOT / 'src/lib/stage_execute.py'),
                   '--stdout', str(self.folder / 'output'),
                   '--stderr', str(self.folder / 'diagnostic'),
                   '--timing', str(self.folder / 'timing'), '--time-command', str(timer),
                   '--', sys.executable, '-c',
                   'import sys; print("stage output"); print("stage diagnostic", file=sys.stderr); sys.exit({})'.format(status)]
        return subprocess.run(command, capture_output=True, text=True)

    def assert_stage(self, timer, status=0, rss='unavailable'):
        result = self.run_stage(timer, status)
        self.assertEqual(result.returncode, status, result.stderr)
        self.assertEqual((self.folder / 'output').read_text(), 'stage output\n')
        self.assertEqual((self.folder / 'diagnostic').read_text(), 'stage diagnostic\n')
        timing = dict(line.split('=', 1) for line in (self.folder / 'timing').read_text().splitlines())
        self.assertGreater(float(timing['elapsed_seconds']), 0)
        self.assertEqual(timing['exit_status'], str(status))
        if rss == 'unavailable':
            self.assertEqual(timing['max_rss_kb'], rss)
        else:
            self.assertGreater(int(timing['max_rss_kb']), 0)

    def test_missing_timer_preserves_stdout_stderr_and_success(self):
        self.assert_stage(self.folder / 'no-such-time')

    def test_non_gnu_timer_is_ignored_and_child_failure_is_preserved(self):
        timer = self.folder / 'bsd-time'
        timer.write_text('#!/bin/sh\n[ "$1" = --version ] && echo "BSD time" && exit 0\nexit 91\n')
        timer.chmod(0o755)
        self.assert_stage(timer, status=17)

    @unittest.skipUnless(Path('/usr/bin/time').exists(), 'GNU time is optional')
    def test_gnu_timer_records_rss(self):
        version = subprocess.run(['/usr/bin/time', '--version'], capture_output=True, text=True)
        if version.returncode or 'GNU Time' not in version.stdout:
            self.skipTest('GNU time is optional')
        self.assert_stage('/usr/bin/time', rss='numeric')


if __name__ == '__main__':
    unittest.main()
