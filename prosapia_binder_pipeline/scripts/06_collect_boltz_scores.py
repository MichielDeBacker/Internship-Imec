#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from prosapia_binder_pipeline.common import (  # noqa: E402
    load_config,
    pair_value,
    read_tsv,
    work_path,
    write_csv,
)

MODEL_RE = re.compile(r"_model_(\d+)\.json$")


def opt_float(data, key):
    value = data.get(key, "")
    if value is None:
        return ""
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    inputs = read_tsv(work_path(cfg, "03_boltz", "boltz_inputs.tsv"))
    mode = str(cfg["scoring"].get("pair_summary_mode", "max"))
    expected_samples = int(cfg["boltz"].get("diffusion_samples", 1))

    all_rows = []
    missing = []

    for item in inputs:
        sequence_id = item["sequence_id"]
        pred_dir = work_path(cfg, "03_boltz", "runs", sequence_id, "predictions", sequence_id)
        confidence_files = sorted(
            glob.glob(str(pred_dir / f"confidence_{sequence_id}_model_*.json"))
        )
        if len(confidence_files) < expected_samples:
            missing.append((sequence_id, len(confidence_files)))
            continue

        for conf_path in confidence_files:
            path = Path(conf_path)
            data = json.loads(path.read_text())
            matrix = data.get("pair_chains_iptm")
            if not isinstance(matrix, dict):
                raise RuntimeError(f"Missing pair_chains_iptm in {path}")

            # Sequence order in 04_prepare_boltz_inputs.py is A, B, C -> indices 0, 1, 2.
            ab, ab_01, ab_10 = pair_value(matrix, 0, 1, mode)
            ac, ac_02, ac_20 = pair_value(matrix, 0, 2, mode)
            bc, bc_12, bc_21 = pair_value(matrix, 1, 2, mode)

            binder_mean = (ac + bc) / 2.0
            binder_min = min(ac, bc)
            binder_delta = abs(ac - bc)

            match = MODEL_RE.search(path.name)
            model_index = int(match.group(1)) if match else -1
            structure = pred_dir / f"{sequence_id}_model_{model_index}.cif"
            if not structure.is_file():
                alt = pred_dir / f"{sequence_id}_model_{model_index}.pdb"
                structure = alt if alt.is_file() else structure

            all_rows.append(
                {
                    "name": sequence_id,
                    "role": item["role"],
                    "design_id": item["design_id"],
                    "backbone_index": item["backbone_index"],
                    "sample": item["sample"],
                    "model_index": model_index,
                    "mpnn_score": item["mpnn_score"],
                    "iptm": opt_float(data, "iptm"),
                    "protein_iptm": opt_float(data, "protein_iptm"),
                    "ptm": opt_float(data, "ptm"),
                    "AB_pair_iptm": ab,
                    "AC_pair_iptm": ac,
                    "BC_pair_iptm": bc,
                    "binder_mean": binder_mean,
                    "binder_min": binder_min,
                    "binder_delta": binder_delta,
                    "AB_A_to_B": ab_01,
                    "AB_B_to_A": ab_10,
                    "AC_A_to_C": ac_02,
                    "AC_C_to_A": ac_20,
                    "BC_B_to_C": bc_12,
                    "BC_C_to_B": bc_21,
                    "confidence_score": opt_float(data, "confidence_score"),
                    "complex_plddt": opt_float(data, "complex_plddt"),
                    "complex_iplddt": opt_float(data, "complex_iplddt"),
                    "complex_pde": opt_float(data, "complex_pde"),
                    "complex_ipde": opt_float(data, "complex_ipde"),
                    "binder_sequence": item["binder_sequence"],
                    "confidence_json": str(path),
                    "predicted_structure": str(structure),
                    "pair_summary_mode": mode,
                }
            )

    if missing:
        preview = "\n".join(f"{name}: {n}" for name, n in missing[:20])
        raise RuntimeError(
            f"Missing/incomplete Boltz outputs for {len(missing)} sequence inputs. "
            f"First entries:\n{preview}"
        )

    if not all_rows:
        raise RuntimeError("No Boltz scores were collected")

    fields = list(all_rows[0].keys())
    all_out = work_path(cfg, "04_results", "boltz_scores_all_models.csv")
    write_csv(all_out, all_rows, fields)

    # One final row per sequence. Boltz orders samples by confidence, but we explicitly
    # select the highest confidence_score so this remains robust to file ordering.
    best_by_name = {}
    for row in all_rows:
        current = best_by_name.get(row["name"])
        score = float(row["confidence_score"]) if row["confidence_score"] != "" else float("-inf")
        cur_score = (
            float(current["confidence_score"])
            if current is not None and current["confidence_score"] != ""
            else float("-inf")
        )
        if current is None or score > cur_score:
            best_by_name[row["name"]] = row

    best_rows = [best_by_name[item["sequence_id"]] for item in inputs]
    best_out = work_path(cfg, "04_results", "boltz_scores_best.csv")
    write_csv(best_out, best_rows, fields)

    compact_fields = [
        "name",
        "role",
        "iptm",
        "AB_pair_iptm",
        "AC_pair_iptm",
        "BC_pair_iptm",
        "binder_mean",
        "binder_min",
        "binder_delta",
    ]
    compact_out = work_path(cfg, "04_results", "cofold_summary.csv")
    compact_rows = []
    numeric_compact = set(compact_fields) - {"name", "role"}
    for row in best_rows:
        compact = {}
        for key in compact_fields:
            value = row[key]
            if key in numeric_compact and value != "":
                value = f"{float(value):.4f}"
            compact[key] = value
        compact_rows.append(compact)
    write_csv(compact_out, compact_rows, compact_fields)

    print(f"Pair summary mode: {mode}")
    print(f"All model rows: {len(all_rows)} -> {all_out}")
    print(f"Best rows: {len(best_rows)} -> {best_out}")
    print(f"Prosapia-style compact table: {compact_out}")


if __name__ == "__main__":
    main()
