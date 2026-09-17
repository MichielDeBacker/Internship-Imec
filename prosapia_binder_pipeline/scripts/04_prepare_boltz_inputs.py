#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from prosapia_binder_pipeline.common import load_config, read_tsv, work_path, write_tsv  # noqa: E402


def add_msa(protein: dict, value):
    if value is None:
        return
    value = str(value)
    if value.lower() == "server":
        return
    if value.lower() == "empty":
        protein["msa"] = "empty"
    else:
        path = Path(value).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Configured MSA does not exist: {path}")
        protein["msa"] = str(path)


def select_rows(rows, per_backbone):
    if str(per_backbone).lower() == "all":
        return rows

    n = int(per_backbone)
    if n <= 0:
        raise ValueError("boltz.sequences_per_backbone must be 'all' or a positive integer")

    grouped = defaultdict(list)
    for row in rows:
        grouped[row["design_id"]].append(row)

    selected = []
    for design_id in sorted(grouped):
        group = grouped[design_id]

        def score_key(row):
            try:
                return float(row["mpnn_score"])
            except (TypeError, ValueError):
                return math.inf

        group.sort(key=lambda r: (score_key(r), int(r["sample"])))
        selected.extend(group[:n])
    return selected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    rows = read_tsv(work_path(cfg, "02_sequence", "sequences.tsv"))
    if not rows:
        raise RuntimeError("Sequence manifest is empty")

    boltz_cfg = cfg["boltz"]
    selected = select_rows(rows, boltz_cfg.get("sequences_per_backbone", "all"))

    # The user requested all 150 backbones to remain represented.
    represented = {row["design_id"] for row in selected}
    expected_backbones = int(cfg["project"]["expected_backbones"])
    if len(represented) != expected_backbones:
        raise RuntimeError(
            f"Boltz selection contains {len(represented)} unique backbones; "
            f"expected {expected_backbones}."
        )

    input_dir = work_path(cfg, "03_boltz", "inputs")
    input_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    for boltz_index, row in enumerate(selected):
        protein_a = {
            "id": "A",
            "sequence": row["chain_A_sequence"],
        }
        protein_b = {
            "id": "B",
            "sequence": row["chain_B_sequence"],
        }
        protein_c = {
            "id": "C",
            "sequence": row["binder_sequence"],
        }
        add_msa(protein_a, boltz_cfg.get("target_a_msa", "empty"))
        add_msa(protein_b, boltz_cfg.get("target_b_msa", "empty"))
        add_msa(protein_c, boltz_cfg.get("binder_msa", "empty"))

        template = {
            "pdb": row["target_template_pdb"],
            "chain_id": ["A", "B"],
            # Boltz PDB templates use subchain ids such as A1 and B1.
            "template_id": ["A1", "B1"],
        }
        if bool(boltz_cfg.get("template_force", False)):
            template["force"] = True
        if boltz_cfg.get("template_threshold_angstrom") is not None:
            template["threshold"] = float(boltz_cfg["template_threshold_angstrom"])

        payload = {
            "version": 1,
            "sequences": [
                {"protein": protein_a},
                {"protein": protein_b},
                {"protein": protein_c},
            ],
            "templates": [template],
        }

        yaml_path = input_dir / f"{row['sequence_id']}.yaml"
        with yaml_path.open("w") as handle:
            yaml.safe_dump(payload, handle, sort_keys=False)

        manifest_rows.append(
            {
                "boltz_index": boltz_index,
                "sequence_index": row["sequence_index"],
                "sequence_id": row["sequence_id"],
                "design_id": row["design_id"],
                "backbone_index": row["backbone_index"],
                "role": row["role"],
                "sample": row["sample"],
                "mpnn_score": row["mpnn_score"],
                "binder_sequence": row["binder_sequence"],
                "boltz_yaml": str(yaml_path),
                "target_template_pdb": row["target_template_pdb"],
            }
        )

    manifest = work_path(cfg, "03_boltz", "boltz_inputs.tsv")
    write_tsv(manifest, manifest_rows, list(manifest_rows[0].keys()))

    uses_server = any(
        str(boltz_cfg.get(key, "empty")).lower() == "server"
        for key in ("target_a_msa", "target_b_msa", "binder_msa")
    )

    print(f"Unique backbones represented: {len(represented)}")
    print(f"Boltz inputs: {len(manifest_rows)}")
    print(f"Uses MSA server: {uses_server}")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
