import json
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class PipelineSmokeTests(unittest.TestCase):

    def test_project_config(self):
        p = ROOT / "pipeline" / "project.json"

        with p.open() as f:
            cfg = json.load(f)

        self.assertIn("design2", cfg["designs"])
        self.assertIn("d2_s0", cfg["primary_candidates"])
        self.assertEqual(
            cfg["rfd3"]["authorized"],
            False,
        )

    def test_cli_help(self):
        p = subprocess.run(
            [
                sys.executable,
                str(ROOT / "pipeline" / "run.py"),
                "--help",
            ],
            text=True,
            capture_output=True,
        )

        self.assertEqual(p.returncode, 0)
        self.assertIn("bootstrap", p.stdout)
        self.assertIn("rfd3-plan", p.stdout)


if __name__ == "__main__":
    unittest.main()
