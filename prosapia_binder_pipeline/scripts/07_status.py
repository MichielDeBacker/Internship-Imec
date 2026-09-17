#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from prosapia_binder_pipeline.common import load_config, read_tsv, work_path  # noqa: E402


def count_rows(path: Path) -> int:
    return len(read_tsv(path)) if path.is_file() else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    cfg = load_config(args.config)

    backbones = work_path(cfg, "01_backbones", "backbones.tsv")
    sequences = work_path(cfg, "02_sequence", "sequences.tsv")
    boltz_inputs = work_path(cfg, "03_boltz", "boltz_inputs.tsv")
    results = work_path(cfg, "04_results", "cofold_summary.csv")

    print(f"backbones.tsv:    {count_rows(backbones)} rows")
    print(f"sequences.tsv:    {count_rows(sequences)} rows")
    print(f"boltz_inputs.tsv: {count_rows(boltz_inputs)} rows")

    seq_done = glob.glob(str(work_path(cfg, "02_sequence", "runs", "*", "DONE.json")))
    boltz_done = glob.glob(str(work_path(cfg, "03_boltz", "runs", "*", "DONE.json")))
    print(f"sequence DONE:    {len(seq_done)}")
    print(f"boltz DONE:       {len(boltz_done)}")
    print(f"cofold summary:   {'present' if results.is_file() else 'not collected yet'}")


if __name__ == "__main__":
    main()
