import json
from pathlib import Path
import shutil
import subprocess
import unittest

REPO = Path(__file__).resolve().parents[1]

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


if __name__ == "__main__":
    unittest.main()
