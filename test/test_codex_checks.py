import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch


REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("codex_checks", REPO / "scripts/codex_checks.py")
checks = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checks)


class CodexChecksTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="codex checks ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name).resolve()
        (self.root / "src/count").mkdir(parents=True)
        (self.root / "test").mkdir()
        shutil.copytree(REPO / "test/fixtures/pipeline_counts",
                        self.root / "test/fixtures/pipeline_counts")

    def event(self, patch_text):
        return {"hook_event_name": "PostToolUse", "tool_name": "apply_patch",
                "cwd": str(self.root), "tool_input": {"command": patch_text}}

    def test_multiple_changes_move_delete_and_spaces(self):
        event = self.event("*** Begin Patch\n*** Add File: one space.sh\n+echo ok\n"
                           "*** Update File: src/old.py\n*** Move to: src/count/new.py\n"
                           "*** Delete File: gone.sh\n*** End Patch")
        self.assertEqual(checks.changed_paths(event, self.root), [
            Path("gone.sh"), Path("one space.sh"), Path("src/count/new.py"), Path("src/old.py")])
        self.assertTrue(checks.needs_units(checks.changed_paths(event, self.root)))

    def test_absolute_subdirectory_and_external_paths(self):
        event = self.event("*** Update File: count/a.py\n*** Add File: ../../outside.sh")
        event["cwd"] = str(self.root / "src")
        self.assertEqual(checks.changed_paths(event, self.root), [Path("src/count/a.py")])
        event = self.event("*** Update File: " + str(self.root / "a.sh"))
        self.assertEqual(checks.changed_paths(event, self.root), [Path("a.sh")])
        (self.root / "escape.sh").symlink_to("/tmp/external.sh")
        self.assertEqual(checks.changed_paths(self.event("*** Update File: escape.sh"), self.root), [])

    def test_trigger_matrix(self):
        for name, units, pipeline in [
            ("README.md", False, False), ("test/new.py", False, False),
            ("scripts/example.sh", False, False), ("config/a.yaml", False, True),
            ("src/deep/a.sh", False, True), ("src/count/a.py", True, True),
            ("src/statistics/a.py", True, False), ("src/metrics/a.py", True, False),
            ("src/visualize/a.R", False, True), ("src/align/a.txt", False, True),
        ]:
            with self.subTest(name=name):
                paths = [Path(name)]
                self.assertEqual(checks.needs_units(paths), units)
                self.assertEqual(checks.needs_pipeline(paths), pipeline)

    def test_shell_syntax_success_and_failure(self):
        path = self.root / "has space.sh"
        event = self.event("*** Add File: has space.sh")
        path.write_text("echo ok\n")
        self.assertIn("PASSED", checks.check_edit(event, self.root)["hookSpecificOutput"]["additionalContext"])
        path.write_text("if then\n")
        result = checks.check_edit(event, self.root)
        self.assertEqual(result["decision"], "block")
        self.assertIn("syntax error", result["reason"])

    def test_deleted_shell_skips_syntax_but_schedules_pipeline(self):
        with patch.object(checks, "run_command") as run, \
                patch.object(checks, "enqueue_pipeline", return_value="queued") as queue:
            checks.check_edit(self.event("*** Delete File: src/gone.sh"), self.root)
        run.assert_not_called()
        queue.assert_called_once_with(self.root)

    def test_units_run_once_failure_keeps_pipeline_trigger(self):
        def fail(command, cwd, log, timeout):
            self.assertEqual(cwd, self.root / "test")
            self.assertEqual(command[-3:], list(checks.UNIT_TESTS))
            log.write_text("specific assertion failure\n")
            return 1
        with patch.object(checks, "run_command", side_effect=fail) as run, \
                patch.object(checks, "enqueue_pipeline", return_value="queued") as queue:
            result = checks.check_edit(self.event("*** Update File: src/count/a.py\n"
                                                 "*** Update File: src/count/b.py"), self.root)
        run.assert_called_once()
        queue.assert_called_once()
        self.assertEqual(result["decision"], "block")
        self.assertIn("specific assertion failure", result["reason"])

    def test_irrelevant_tool_and_file(self):
        event = self.event("*** Add File: a.sh")
        event["tool_name"] = "Bash"
        self.assertEqual(checks.check_edit(event, self.root), {})
        self.assertEqual(checks.check_edit(self.event("*** Update File: README.md"), self.root), {})

    def test_timeout_and_missing_executable(self):
        log = self.root / "timeout.log"
        self.assertEqual(checks.run_command([sys.executable, "-c", "import time; time.sleep(20)"],
                                           self.root, log, 0.05), 124)
        self.assertIn("TIMEOUT", log.read_text())
        self.assertEqual(checks.run_command([str(self.root / "missing")], self.root, log, 1), 127)

    def wait_for(self, condition):
        deadline = time.monotonic() + 15
        while time.monotonic() < deadline:
            if condition():
                return
            time.sleep(0.03)
        self.fail("Worker did not reach expected state: " + repr(checks.read_state(self.root)))

    def test_background_worker_coalesces_and_preserves_existing_genomes(self):
        existing = self.root / "test_genomes"
        existing.mkdir()
        (existing / "keep").write_text("existing data")
        # A real child process exercises inherited flock descriptors through RTK.
        (self.root / "src/sbst_fromDwl.sh").write_text(
            "#!/bin/bash\npython3 - \"$@\" <<'PY'\n"
            "from pathlib import Path\nimport time\n"
            "import shutil\nimport sys\n"
            "date = sys.argv[1]\n"
            "results = Path(sys.argv[sys.argv.index('--out-dir') + 1])\n"
            "statistics = results / 'Oikdio1_Oikalb2_Oikvan3' / date / 'statistics'\n"
            "shutil.copytree('test/fixtures/pipeline_counts', statistics)\n"
            "for path in statistics.glob('*/*/*.tsv'):\n"
            "    path.rename(path.with_name(path.name.replace('20260420', date)))\n"
            "p = Path('runs')\nn = int(p.read_text()) + 1 if p.exists() else 1\n"
            "p.write_text(str(n))\n"
            "if n == 1:\n"
            "    deadline = time.monotonic() + 10\n"
            "    while not Path('release').exists() and time.monotonic() < deadline:\n"
            "        time.sleep(0.02)\nPY\n")
        checks.enqueue_pipeline(self.root)
        self.wait_for(lambda: (self.root / "runs").exists())
        checks.enqueue_pipeline(self.root)
        checks.enqueue_pipeline(self.root)
        self.assertEqual((self.root / "runs").read_text(), "1")
        (self.root / "release").touch()
        self.wait_for(lambda: checks.read_state(self.root).get("status") == "passed")
        self.assertEqual((self.root / "runs").read_text(), "2")
        state = checks.read_state(self.root)
        self.assertEqual(state["generation"], 3)
        self.assertEqual(state["requested"], 3)
        self.assertEqual(state["counting_comparison"]["matched_count"], 8)
        self.assertEqual((existing / "keep").read_text(), "existing data")
        logs = list(checks.state_dir(self.root).glob("pipeline-*/pipeline.log"))
        self.assertEqual(len(logs), 2)
        for log in logs:
            self.assertIn("EXIT CODE: 0", log.read_text())
            self.assertIn(str(log.parent / "genomes"), log.read_text())
            self.assertIn("Counting TSVs: 8/8 match", log.read_text())
            result = json.loads((log.parent / "result.json").read_text())
            self.assertEqual(result["counting_comparison"]["status"], "passed")

    def test_worker_rejects_successful_pipeline_without_tsvs(self):
        (self.root / "src/sbst_fromDwl.sh").write_text("#!/bin/bash\nexit 0\n")
        checks.enqueue_pipeline(self.root)
        self.wait_for(lambda: checks.read_state(self.root).get("status") == "failed")
        state = checks.read_state(self.root)
        self.assertEqual(state["pipeline_exit_code"], 0)
        self.assertEqual(state["exit_code"], 1)
        self.assertEqual(state["counting_comparison"]["matched_count"], 0)
        self.assertEqual(len(state["counting_comparison"]["files"]), 8)

    def test_worker_timeout_skips_comparison_and_preserves_exit_code(self):
        with checks.state_lock(self.root):
            checks.write_state(self.root, dict(requested=1))
        lock_path = checks.state_dir(self.root) / "worker.lock"
        # pipeline_worker owns and closes its descriptor.
        lock_fd = os.open(lock_path, os.O_CREAT | os.O_WRONLY, 0o600)

        def timeout(command, cwd, log, limit):
            self.assertEqual(limit, 600)
            log.write_text("TIMEOUT\n")
            return 124

        with patch.object(checks, "run_command", side_effect=timeout), \
                patch.object(checks, "compare_counting_tsvs") as compare:
            checks.pipeline_worker(self.root, lock_fd)
        compare.assert_not_called()
        state = checks.read_state(self.root)
        self.assertEqual(state["status"], "failed")
        self.assertEqual(state["exit_code"], 124)
        self.assertEqual(state["counting_comparison"]["status"], "skipped")

    def test_worker_records_tsv_mismatch_as_failure_after_pipeline_success(self):
        checks.write_state(self.root, dict(requested=1))
        lock_fd = os.open(checks.state_dir(self.root) / "worker.lock",
                          os.O_CREAT | os.O_WRONLY, 0o600)

        def pipeline(command, cwd, log, limit):
            date = command[4]
            results = Path(command[command.index("--out-dir") + 1])
            statistics = results / checks.COUNT_TRIO / date / "statistics"
            shutil.copytree(self.root / "test/fixtures/pipeline_counts", statistics)
            for path in statistics.glob("*/*/*.tsv"):
                path.rename(path.with_name(path.name.replace("20260420", date)))
            actual = sorted(statistics.glob("*/*/*.tsv"))[0]
            actual.write_bytes(actual.read_bytes() + b"\n")
            log.write_text("EXIT CODE: 0\n")
            return 0

        with patch.object(checks, "run_command", side_effect=pipeline):
            checks.pipeline_worker(self.root, lock_fd)
        state = checks.read_state(self.root)
        self.assertEqual(state["status"], "failed")
        self.assertEqual(state["pipeline_exit_code"], 0)
        self.assertEqual(state["exit_code"], 1)
        self.assertEqual(state["counting_comparison"]["matched_count"], 7)
        recorded = json.loads((Path(state["run_dir"]) / "result.json").read_text())
        self.assertEqual(recorded["counting_comparison"], state["counting_comparison"])
        self.assertIn("@@", Path(state["log"]).read_text())
        self.assertIn("PIPELINE CHECK: FAILED", Path(state["log"]).read_text())

    def test_worker_failure_recorded(self):
        (self.root / "src/sbst_fromDwl.sh").write_text("#!/bin/bash\nexit 7\n")
        checks.enqueue_pipeline(self.root)
        self.wait_for(lambda: checks.read_state(self.root).get("status") == "failed")
        self.assertEqual(checks.read_state(self.root)["exit_code"], 7)
        self.assertEqual(checks.read_state(self.root)["counting_comparison"]["status"], "skipped")

    def test_spawn_failure_does_not_keep_lock(self):
        with patch.object(checks.subprocess, "Popen", side_effect=OSError("spawn failed")):
            with self.assertRaises(OSError):
                checks.enqueue_pipeline(self.root)
        self.assertEqual(checks.read_state(self.root)["status"], "worker_start_failed")
        self.test_worker_failure_recorded()

    def test_hook_json_and_malformed_event(self):
        command = [sys.executable, str(REPO / "scripts/codex_checks.py"), "--root", str(self.root)]
        result = subprocess.run(command, input=json.dumps(self.event("*** Update File: README.md")),
                                text=True, capture_output=True)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(json.loads(result.stdout), {})
        result = subprocess.run(command, input="not json", text=True, capture_output=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("Codex checks failed", result.stderr)

    def test_configured_command_accepts_native_patch_event(self):
        configuration = json.loads((REPO / ".codex/hooks.json").read_text())
        hook = configuration["hooks"]["PostToolUse"][0]["hooks"][0]
        # Exercise the exact configured shell command, including repository discovery.
        (REPO / "tmp").mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(prefix="hook event ", dir=str(REPO / "tmp")) as folder:
            path = Path(folder) / "sample.sh"
            path.write_text("echo valid\n")
            event = {"tool_name": "apply_patch", "cwd": str(REPO),
                     "tool_input": {"command": "*** Begin Patch\n*** Add File: " +
                                    str(path) + "\n+echo valid\n*** End Patch"}}
            result = subprocess.run(["rtk", "proxy", "bash", "-c", hook["command"]],
                                    input=json.dumps(event), cwd=str(REPO / "test"),
                                    text=True, capture_output=True, check=True)
            context = json.loads(result.stdout)["hookSpecificOutput"]["additionalContext"]
            self.assertIn("PASSED", context)
            self.assertNotIn("Pipeline", context)


@unittest.skipUnless(shutil.which("codex"), "Codex CLI is required for rule validation")
class CodexRulesTests(unittest.TestCase):
    def test_forbidden_and_unrelated_commands(self):
        rules = REPO / ".codex/rules/safety.rules"
        forbidden = [["rm", "-rf", "/unused"], ["git", "push", "--force"],
                     ["git", "push", "-f"], ["git", "reset", "--hard"],
                     ["git", "clean", "-f"], ["git", "clean", "-fd"]]
        for wrapper in ([], ["rtk"], ["rtk", "proxy"]):
            for command in forbidden + [["git", "status"], ["git", "diff"]]:
                with self.subTest(command=wrapper + command):
                    result = subprocess.run(
                        ["rtk", "proxy", "codex", "execpolicy", "check", "--rules", str(rules),
                         "--", *wrapper, *command], capture_output=True, text=True, check=True)
                    self.assertEqual(json.loads(result.stdout).get("decision") == "forbidden",
                                     command in forbidden)


class CountingComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="counting comparison ")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.baseline = self.root / "test/fixtures/pipeline_counts"
        shutil.copytree(REPO / "test/fixtures/pipeline_counts", self.baseline)
        self.results = self.root / "results"
        self.date = "20990102030405"
        statistics = self.results / checks.COUNT_TRIO / self.date / "statistics"
        shutil.copytree(self.baseline, statistics)
        self.actuals = []
        for path in sorted(statistics.glob("*/*/*.tsv")):
            actual = path.with_name(path.name.replace("20260420", self.date))
            path.rename(actual)
            self.actuals.append(actual)
        self.log = self.root / "pipeline.log"

    def compare(self):
        return checks.compare_counting_tsvs(self.root, self.results, self.date, self.log)

    def test_all_eight_historical_files_match_with_different_run_date(self):
        self.assertEqual(len(self.actuals), 8)
        result = self.compare()
        self.assertEqual(result["status"], "passed")
        self.assertEqual(result["matched_count"], 8)
        self.assertEqual(result["expected_count"], 8)
        self.assertIn("Counting TSVs: 8/8 match", self.log.read_text())

    def test_one_changed_count_in_each_file_fails_with_diff(self):
        for actual in self.actuals:
            with self.subTest(file=actual.name):
                original = actual.read_bytes()
                lines = original.splitlines(keepends=True)
                columns = lines[1].split(b"\t")
                columns[1] = str(int(columns[1]) + 1).encode()
                lines[1] = b"\t".join(columns)
                actual.write_bytes(b"".join(lines))
                result = self.compare()
                self.assertEqual(result["status"], "failed")
                self.assertEqual(result["matched_count"], 7)
                failed = [f for f in result["files"] if f["status"] == "failed"]
                self.assertEqual(failed[0]["actual"], str(actual))
                self.assertIn("TSV bytes differ", failed[0]["error"])
                self.assertIn("@@", self.log.read_text())
                actual.write_bytes(original)

    def test_byte_level_differences_are_not_normalized(self):
        actual = self.actuals[0]
        original = actual.read_bytes()
        lines = original.splitlines(keepends=True)
        for label, changed in [
            ("row order", b"".join([lines[0], lines[2], lines[1], *lines[3:]])),
            ("CRLF", original.replace(b"\n", b"\r\n")),
            ("missing final newline", original.rstrip(b"\n")),
            ("changed header", original.replace(b"mutType", b"mut_type", 1)),
        ]:
            with self.subTest(change=label):
                actual.write_bytes(changed)
                self.assertEqual(self.compare()["matched_count"], 7)

    def test_missing_or_empty_files_cannot_pass(self):
        actual = self.actuals[0]
        actual.unlink()
        result = self.compare()
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["matched_count"], 7)
        actual.touch()
        failed = [f for f in self.compare()["files"] if f["status"] == "failed"]
        self.assertIn("Empty", failed[0]["error"])

    def test_missing_or_empty_baseline_cannot_pass(self):
        expected = sorted(self.baseline.glob("*/*/*.tsv"))[0]
        expected.unlink()
        self.assertEqual(self.compare()["matched_count"], 7)
        expected.touch()
        failed = [f for f in self.compare()["files"] if f["status"] == "failed"]
        self.assertIn("Empty", failed[0]["error"])

    def test_non_cds_outputs_are_required(self):
        for actual in self.actuals:
            if actual.name.endswith("_ncds.tsv"):
                actual.unlink()
        result = self.compare()
        self.assertEqual(result["status"], "failed")
        self.assertEqual(result["matched_count"], 4)
        self.assertEqual(result["expected_count"], 8)


if __name__ == "__main__":
    unittest.main()
