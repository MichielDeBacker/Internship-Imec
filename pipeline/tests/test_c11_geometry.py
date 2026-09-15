import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


class C11GeometryTests(unittest.TestCase):

    def test_cli_help(self):
        p = subprocess.run(
            [
                sys.executable,
                str(
                    ROOT
                    / "pipeline"
                    / "lib"
                    / "c11_geometry.py"
                ),
                "--help",
            ],
            text=True,
            capture_output=True,
        )

        self.assertEqual(
            p.returncode,
            0,
        )

        self.assertIn(
            "--parent",
            p.stdout,
        )


if __name__ == "__main__":
    unittest.main()
