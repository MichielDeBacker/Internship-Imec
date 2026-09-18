#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import os
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from prosapia_binder_pipeline.common import (  # noqa: E402
    atomic_json,
    ensure_file,
    expand_path,
    get_indexed_record,
    load_config,
    run,
    work_path,
)


def uses_msa_server(cfg) -> bool:
    boltz = cfg["boltz"]
    return any(
        str(boltz.get(key, "empty")).lower() == "server"
        for key in ("target_a_msa", "target_b_msa", "binder_msa")
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--index", required=True, type=int)
    ap.add_argument("--override", action="store_true")
    args = ap.parse_args()

    cfg = load_config(args.config)
    manifest = work_path(cfg, "03_boltz", "boltz_inputs.tsv")
    row = get_indexed_record(manifest, args.index)

    sequence_id = row["sequence_id"]
    input_yaml = ensure_file(row["boltz_yaml"], "Boltz input YAML")
    out_dir = work_path(cfg, "03_boltz", "runs", sequence_id)
    out_dir.mkdir(parents=True, exist_ok=True)

    expected_dir = out_dir / f"boltz_results_{sequence_id}" / "predictions" / sequence_id
    existing = sorted(glob.glob(str(expected_dir / f"confidence_{sequence_id}_model_*.json")))
    expected_samples = int(cfg["boltz"]["diffusion_samples"])
    if len(existing) >= expected_samples and not args.override:
        print(f"Already complete: {sequence_id} ({len(existing)} confidence JSON files)")
        return

    boltz_cfg = cfg["boltz"]
    cmd = [
        str(cfg["external"].get("boltz_executable", "boltz")),
        "predict",
        str(input_yaml),
        "--out_dir", str(out_dir),
        "--recycling_steps", str(boltz_cfg.get("recycling_steps", 3)),
        "--sampling_steps", str(boltz_cfg.get("sampling_steps", 200)),
        "--diffusion_samples", str(expected_samples),
        "--devices", "1",
        "--accelerator", "gpu",
    ]

    if args.override:
        cmd.append("--override")
    if bool(boltz_cfg.get("use_potentials", False)):
        cmd.append("--use_potentials")
    if bool(boltz_cfg.get("no_kernels", False)):
        cmd.append("--no_kernels")
    if bool(boltz_cfg.get("write_full_pae", False)):
        cmd.append("--write_full_pae")
    if bool(boltz_cfg.get("write_full_pde", False)):
        cmd.append("--write_full_pde")
    if uses_msa_server(cfg):
        cmd.extend(["--use_msa_server", "--msa_server_url", str(boltz_cfg["msa_server_url"])])

    env = os.environ.copy()
    cache = expand_path(cfg["external"].get("boltz_cache", work_path(cfg, "cache", "boltz")))
    cache.mkdir(parents=True, exist_ok=True)
    env["BOLTZ_CACHE"] = str(cache)

    run(cmd, env=env)

    confidence = sorted(glob.glob(str(expected_dir / f"confidence_{sequence_id}_model_*.json")))
    if len(confidence) < expected_samples:
        raise RuntimeError(
            f"{sequence_id}: expected {expected_samples} confidence JSON files, found {len(confidence)}"
        )

    atomic_json(
        out_dir / "DONE.json",
        {
            "stage": "boltz",
            "boltz_index": args.index,
            "sequence_id": sequence_id,
            "input_yaml": str(input_yaml),
            "output_dir": str(out_dir),
            "confidence_files": confidence,
            "command": shlex.join(cmd),
            "boltz_cache": str(cache),
        },
    )
    print(f"Completed {sequence_id}: {len(confidence)} diffusion sample(s)")


if __name__ == "__main__":
    main()
