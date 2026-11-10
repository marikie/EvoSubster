#!/usr/bin/env python3
"""Execute one stage and report resources even when GNU time is unavailable."""

import argparse
from pathlib import Path
import subprocess
import tempfile
import time


def is_gnu_time(command):
    try:
        version = subprocess.run([command, '--version'], capture_output=True,
                                 text=True, timeout=5)
        return version.returncode == 0 and 'GNU Time' in version.stdout
    except (OSError, subprocess.SubprocessError):
        return False


def execute(command, output, diagnostic, timing, time_command='/usr/bin/time'):
    """Preserve child output and status; timing support never gates execution."""
    max_rss = 'unavailable'
    with tempfile.TemporaryDirectory(prefix='evosubster-time-') as temporary:
        rss_path = Path(temporary) / 'rss'
        argv = command
        if is_gnu_time(time_command):
            argv = [time_command, '-f', '%M', '-o', str(rss_path), '--', *command]
        started = time.monotonic()
        with open(output, 'wb') as stdout, open(diagnostic, 'wb') as stderr:
            try:
                status = subprocess.run(argv, stdout=stdout, stderr=stderr).returncode
                if status < 0:
                    status = 128 - status
            except OSError as error:
                stderr.write(('Cannot execute stage: {}\n'.format(error)).encode())
                status = 127
        elapsed = time.monotonic() - started
        if rss_path.exists():
            # GNU time can prepend a diagnostic when the child exits nonzero.
            values = rss_path.read_text().splitlines()
            if values and values[-1].isdigit():
                max_rss = values[-1]
    Path(timing).write_text('elapsed_seconds={:.6f}\nmax_rss_kb={}\nexit_status={}\n'.format(
        elapsed, max_rss, status))
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stdout', required=True)
    parser.add_argument('--stderr', required=True)
    parser.add_argument('--timing', required=True)
    parser.add_argument('--time-command', default='/usr/bin/time')
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('a stage command is required')
    return execute(command, args.stdout, args.stderr, args.timing, args.time_command)


if __name__ == '__main__':
    raise SystemExit(main())
