#!/usr/bin/env python3
"""Codex PostToolUse checks and a serialized, coalescing pipeline worker."""

import argparse
from contextlib import contextmanager
import difflib
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import threading
import time


ROOT = Path(__file__).resolve().parents[1]
UNIT_TESTS = (
    "test_single_sbst_2TSVs", "test_single_sbst_2TSVs_errprob", "test_subRatio"
)
ACCESSIONS = ("GCA_907165135.1", "GCA_004367875.1", "GCA_004367855.1")
COUNT_BASELINE_DATE = "20260420"
COUNT_TRIO = "Oikdio1_Oikalb2_Oikvan3"


def counting_tsv_paths(date):
    """The eight required counting outputs, relative to statistics/."""
    for accession, organism in zip(ACCESSIONS[1:], ("Oikalb2", "Oikvan3")):
        for kind, suffix in (("singlenuc", ""), ("dinuc", "_dinuc")):
            for region in ("", "_ncds"):
                yield Path(organism) / kind / ("{}_{}_{}{}{}.tsv".format(
                    accession, organism, date, suffix, region))


def compare_counting_tsvs(root, results_dir, date, log):
    """Compare bytes against the saved historical TSVs, without normalization."""
    baseline = root / "test" / "fixtures" / "pipeline_counts"
    statistics = results_dir / COUNT_TRIO / date / "statistics"
    files = []
    with log.open("ab") as output:
        output.write(b"\nCounting TSV comparison (byte-for-byte):\n")
        for expected_relative, actual_relative in zip(
                counting_tsv_paths(COUNT_BASELINE_DATE), counting_tsv_paths(date)):
            expected = baseline / expected_relative
            actual = statistics / actual_relative
            item = dict(expected=str(expected), actual=str(actual))
            try:
                expected_bytes = expected.read_bytes()
                actual_bytes = actual.read_bytes()
            except OSError as error:
                item.update(status="failed", error=str(error))
            else:
                if not expected_bytes or not actual_bytes:
                    item.update(status="failed", error="Empty expected or actual TSV")
                elif expected_bytes == actual_bytes:
                    item["status"] = "passed"
                else:
                    offset = next((i for i, (left, right) in enumerate(
                        zip(expected_bytes, actual_bytes)) if left != right),
                        min(len(expected_bytes), len(actual_bytes)))
                    item.update(status="failed", error=(
                        "TSV bytes differ at offset {}; expected {} bytes, actual {} bytes"
                        .format(offset, len(expected_bytes), len(actual_bytes))))
            output.write(("{}: {}{}\n".format(
                item["status"].upper(), actual,
                ": " + item["error"] if "error" in item else "")).encode())
            if item.get("error", "").startswith("TSV bytes differ"):
                for line in difflib.diff_bytes(
                        difflib.unified_diff, expected_bytes.splitlines(keepends=True),
                        actual_bytes.splitlines(keepends=True),
                        fromfile=os.fsencode(expected), tofile=os.fsencode(actual)):
                    output.write(line)
                    if not line.endswith(b"\n"):
                        output.write(b"\n\\ No newline at end of file\n")
            files.append(item)
        matched = sum(item["status"] == "passed" for item in files)
        output.write("Counting TSVs: {}/{} match\n".format(matched, len(files)).encode())
    return dict(status="passed" if matched == len(files) else "failed",
                expected_count=len(files), matched_count=matched, files=files)


def state_dir(root):
    path = root / "tmp" / "codex-checks"
    path.mkdir(parents=True, exist_ok=True)
    return path


@contextmanager
def state_lock(root):
    with (state_dir(root) / "state.lock").open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def read_state(root):
    path = state_dir(root) / "status.json"
    return json.loads(path.read_text()) if path.exists() else {}


def write_state(root, state):
    path = state_dir(root) / "status.json"
    temporary = path.with_suffix(".new")
    temporary.write_text(json.dumps(state, indent=2) + "\n")
    temporary.replace(path)


def changed_paths(event, root):
    """Read native Codex patches; retain deleted/source paths for test triggers."""
    if event.get("tool_name") not in ("apply_patch", "Edit", "Write"):
        return []
    data = event.get("tool_input", {})
    if not isinstance(data, dict):
        raise ValueError("Expected an object in tool_input")
    names = []
    if data.get("file_path"):
        names.append(data["file_path"])
    patch = data.get("command", data.get("patch", ""))
    for line in patch.splitlines():
        for prefix in ("*** Add File: ", "*** Update File: ",
                       "*** Delete File: ", "*** Move to: "):
            if line.startswith(prefix):
                names.append(line[len(prefix):])
                break
    cwd = Path(event.get("cwd", str(root)))
    paths = set()
    for name in names:
        path = Path(name)
        path = (cwd / path).resolve() if not path.is_absolute() else path.resolve()
        try:
            paths.add(path.relative_to(root.resolve()))
        except ValueError:
            continue  # Never run project checks on files outside this repository.
    return sorted(paths)


def needs_units(paths):
    return any(p.suffix == ".py" and p.parts[:2] in
               (("src", "count"), ("src", "statistics"), ("src", "metrics"))
               for p in paths)


def needs_pipeline(paths):
    # Bash's original */src/*.sh pattern also matches nested shell scripts.
    return any(p.parts[0] == "config" or
               (p.parts[0] == "src" and p.suffix == ".sh") or
               p.parts[:2] in (("src", "align"), ("src", "count"),
                              ("src", "visualize")) for p in paths)


def run_command(command, cwd, log, timeout):
    """Terminate the whole child process group on timeout, including LAST jobs."""
    with log.open("w") as output:
        output.write("Command: " + repr(command) + "\n")
        output.flush()
        try:
            process = subprocess.Popen(command, cwd=str(cwd), stdout=output,
                                       stderr=subprocess.STDOUT, start_new_session=True)
        except OSError as error:
            output.write(str(error) + "\n")
            return 127
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
            # Allow cleanup, then kill descendants even if the group leader exited.
            try:
                process.wait(timeout=2)
            except subprocess.TimeoutExpired:
                pass
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            output.write("TIMEOUT\n")
            code = 124
        output.write("EXIT CODE: {}\n".format(code))
        return code


def enqueue_pipeline(root):
    with state_lock(root):
        state = read_state(root)
        state["requested"] = state.get("requested", 0) + 1
        write_state(root, state)
        worker_lock = (state_dir(root) / "worker.lock").open("a")
        try:
            try:
                fcntl.flock(worker_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                return "Pipeline rerun queued; status: " + str(state_dir(root) / "status.json")
            state["status"] = "queued"
            write_state(root, state)
            with (state_dir(root) / "worker.log").open("a") as output:
                process = subprocess.Popen(
                    ["rtk", "proxy", sys.executable, str(Path(__file__).resolve()),
                     "--worker", "--root", str(root), "--lock-fd", str(worker_lock.fileno())],
                    cwd=str(root), stdin=subprocess.DEVNULL, stdout=output, stderr=output,
                    start_new_session=True, pass_fds=(worker_lock.fileno(),))
                # Reap children in long-lived callers without waiting in the hook.
                threading.Thread(target=process.wait, daemon=True).start()
        except Exception:
            state["status"] = "worker_start_failed"
            write_state(root, state)
            raise
        finally:
            # Closing our copy leaves the inherited worker descriptor locked.
            worker_lock.close()
    return "Pipeline queued; status: " + str(state_dir(root) / "status.json")


def pipeline_worker(root, lock_fd):
    try:
        while True:
            with state_lock(root):
                state = read_state(root)
                generation = state["requested"]
                for key in ("exit_code", "pipeline_exit_code", "counting_comparison",
                            "finished_at", "error"):
                    state.pop(key, None)
                run_dir = Path(tempfile.mkdtemp(prefix="pipeline-", dir=str(state_dir(root))))
                log = run_dir / "pipeline.log"
                state.update(status="running", generation=generation, log=str(log),
                             run_dir=str(run_dir), started_at=time.time())
                write_state(root, state)
            date = time.strftime("%Y%m%d%H%M%S")
            command = ["rtk", "proxy", "bash", str(root / "src" / "sbst_fromDwl.sh"),
                       date, *ACCESSIONS,
                       "--genome-dir", str(run_dir / "genomes"),
                       "--out-dir", str(run_dir / "results"),
                       "--train-cache-dir", str(run_dir / "train-cache")]
            code = run_command(command, root, log, 600)
            pipeline_code = code
            if code == 0:
                comparison = compare_counting_tsvs(root, run_dir / "results", date, log)
                if comparison["status"] != "passed":
                    code = 1
            else:
                comparison = dict(status="skipped", reason="Pipeline did not exit successfully")
                with log.open("a") as output:
                    output.write("Counting TSV comparison SKIPPED: pipeline exit {}\n".format(code))
            result = dict(generation=generation, exit_code=code,
                          pipeline_exit_code=pipeline_code, counting_comparison=comparison,
                          log=str(log),
                          status="passed" if code == 0 else "failed", finished_at=time.time())
            with log.open("a") as output:
                output.write("PIPELINE CHECK: {} (exit {})\n".format(result["status"].upper(), code))
            (run_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n")
            with state_lock(root):
                state = read_state(root)
                state.update(result)
                if state["requested"] == generation:
                    write_state(root, state)
                    # Serialize idle transition with enqueue, avoiding lost requests.
                    fcntl.flock(lock_fd, fcntl.LOCK_UN)
                    return
                state["status"] = "queued"
                write_state(root, state)
    except Exception as error:
        with state_lock(root):
            state = read_state(root)
            state.update(status="worker_failed", error=str(error))
            write_state(root, state)
        raise
    finally:
        os.close(lock_fd)


def check_edit(event, root):
    paths = changed_paths(event, root)
    if not paths:
        return {}
    messages = []
    failed = False
    checks = []
    for path in paths:
        if path.suffix == ".sh" and (root / path).is_file():
            checks.append(("bash syntax: " + str(path),
                           ["rtk", "proxy", "bash", "-n", str(root / path)], root, 15))
    if needs_units(paths):
        checks.append(("Python unit tests", ["rtk", "proxy", "python3", "-m", "unittest",
                                            *UNIT_TESTS], root / "test", 120))
    # Serialize unit suites: the legacy tests share result fixture filenames.
    with (state_dir(root) / "checks.lock").open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        for label, command, cwd, timeout in checks:
            fd, name = tempfile.mkstemp(prefix="check-", suffix=".log", dir=str(state_dir(root)))
            os.close(fd)
            log = Path(name)
            code = run_command(command, cwd, log, timeout)
            messages.append("{}: {} (exit {}); log: {}".format(
                label, "PASSED" if code == 0 else "FAILED", code, log))
            if code:
                failed = True
                messages.append(log.read_text(errors="replace")[-6000:])
    if needs_pipeline(paths):
        messages.append(enqueue_pipeline(root))
    if not messages:
        return {}
    result = {"hookSpecificOutput": {"hookEventName": "PostToolUse",
                                     "additionalContext": "\n".join(messages)}}
    if failed:
        result.update(decision="block", reason="\n".join(messages))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--lock-fd", type=int)
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--queue-pipeline", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    if args.worker:
        if args.lock_fd is None:
            parser.error("--worker requires --lock-fd")
        pipeline_worker(root, args.lock_fd)
    elif args.status:
        with state_lock(root):
            print(json.dumps(read_state(root), indent=2))
    elif args.queue_pipeline:
        print(enqueue_pipeline(root))
    else:
        try:
            print(json.dumps(check_edit(json.load(sys.stdin), root)))
        except Exception as error:
            print("Codex checks failed: " + str(error), file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
