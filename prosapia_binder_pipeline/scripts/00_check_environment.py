#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from prosapia_binder_pipeline.common import expand_path, load_config  # noqa: E402


def git_info(path: Path):
    if not (path / ".git").exists():
        return None, None
    commit = subprocess.check_output(
        ["git", "-C", str(path), "rev-parse", "HEAD"], text=True
    ).strip()
    status = subprocess.check_output(
        ["git", "-C", str(path), "status", "--porcelain"], text=True
    ).strip()
    return commit, status


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    failures = []

    prosapia = expand_path(cfg["external"]["prosapia_root"])
    print(f"Prosapia root (READ ONLY by this pipeline): {prosapia}")
    if prosapia.exists():
        commit, status = git_info(prosapia)
        if commit:
            print(f"  git commit: {commit}")
            print(f"  pre-existing dirty state: {'yes' if status else 'no'}")
    else:
        print("  note: checkout not found; pipeline can still run because it does not patch Prosapia")

    mpnn_root = expand_path(cfg["external"]["proteinmpnn_root"])
    mpnn = mpnn_root / "protein_mpnn_run.py"
    print(f"ProteinMPNN: {mpnn}")
    if not mpnn.is_file():
        failures.append(f"Missing ProteinMPNN entrypoint: {mpnn}")

    mpnn_python = str(cfg["external"].get("proteinmpnn_python", "python"))
    if "/" in mpnn_python:
        ok = Path(mpnn_python).is_file()
    else:
        ok = shutil.which(mpnn_python) is not None
    print(f"ProteinMPNN Python: {mpnn_python} ({'OK' if ok else 'NOT FOUND'})")
    if not ok:
        failures.append(f"ProteinMPNN Python not found: {mpnn_python}")

    boltz = str(cfg["external"].get("boltz_executable", "boltz"))
    if "/" in boltz:
        boltz_ok = Path(boltz).is_file()
    else:
        boltz_ok = shutil.which(boltz) is not None
    print(f"Boltz executable: {boltz} ({'OK' if boltz_ok else 'NOT FOUND in current shell'})")
    if not boltz_ok:
        print("  This may still be fine if Boltz is available only inside the Slurm job environment.")

    backbone_glob = cfg["paths"]["backbone_glob"]
    print(f"Backbone glob: {backbone_glob}")
    print(f"Work root: {cfg['paths']['work_root']}")

    if failures:
        print("\nFAILURES:")
        for failure in failures:
            print(f"- {failure}")
        raise SystemExit(2)

    print("\nEnvironment check passed for the required sequencing files.")


if __name__ == "__main__":
    main()
