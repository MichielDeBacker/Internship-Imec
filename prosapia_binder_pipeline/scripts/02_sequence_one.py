#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from prosapia_binder_pipeline.common import (  # noqa: E402
    atomic_json,
    ensure_file,
    expand_path,
    fasta_records,
    get_indexed_record,
    load_config,
    run,
    work_path,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--index", required=True, type=int)
    args = ap.parse_args()

    cfg = load_config(args.config)
    manifest = work_path(cfg, "01_backbones", "backbones.tsv")
    row = get_indexed_record(manifest, args.index)

    design_id = row["design_id"]
    pdb = ensure_file(row["canonical_pdb"], "canonical PDB")

    mpnn_root = expand_path(cfg["external"]["proteinmpnn_root"])
    mpnn_script = ensure_file(mpnn_root / "protein_mpnn_run.py", "ProteinMPNN entrypoint")
    mpnn_python = str(cfg["external"].get("proteinmpnn_python", "python"))

    seq_cfg = cfg["sequencing"]
    nseq = int(seq_cfg["num_sequences_per_backbone"])
    batch_size = int(seq_cfg["batch_size"])
    if nseq % batch_size != 0:
        raise ValueError(
            "ProteinMPNN requires num_sequences_per_backbone to be divisible by batch_size "
            f"with this runner; got {nseq} and {batch_size}."
        )

    out_dir = work_path(cfg, "02_sequence", "runs", design_id)
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = out_dir / "seqs" / f"{pdb.stem}.fa"
    done = out_dir / "DONE.json"

    if fasta.is_file():
        records = fasta_records(fasta)
        sampled = records[1:] if records else []
        if len(sampled) == nseq:
            print(f"Already complete: {design_id} ({len(sampled)} sampled sequences)")
            return

    cmd = [
        mpnn_python,
        str(mpnn_script),
        "--pdb_path", str(pdb),
        "--pdb_path_chains", "C",
        "--out_folder", str(out_dir),
        "--model_name", str(seq_cfg["model_name"]),
        "--num_seq_per_target", str(nseq),
        "--batch_size", str(batch_size),
        "--sampling_temp", str(seq_cfg["sampling_temperature"]),
        "--backbone_noise", str(seq_cfg["backbone_noise"]),
        "--save_score", "1",
        "--seed", str(int(seq_cfg["seed_base"]) + args.index),
    ]
    if bool(seq_cfg.get("use_soluble_model", False)):
        cmd.append("--use_soluble_model")

    run(cmd)

    ensure_file(fasta, "ProteinMPNN FASTA")
    records = fasta_records(fasta)
    sampled = records[1:]
    if len(sampled) != nseq:
        raise RuntimeError(
            f"{design_id}: expected {nseq} sampled sequences after the native record, "
            f"found {len(sampled)} in {fasta}"
        )

    for _, sequence in sampled:
        if "/" in sequence:
            raise RuntimeError(
                f"{design_id}: output contains multiple designed chains, but only C should be designed"
            )
        if len(sequence) != int(row["chain_C_length"]):
            raise RuntimeError(
                f"{design_id}: sampled sequence length {len(sequence)} != binder length {row['chain_C_length']}"
            )

    atomic_json(
        done,
        {
            "stage": "sequence",
            "backbone_index": args.index,
            "design_id": design_id,
            "designed_chain": "C",
            "fixed_chains": ["A", "B"],
            "num_sequences": nseq,
            "command": shlex.join(cmd),
            "fasta": str(fasta),
        },
    )
    print(f"Completed {design_id}: {nseq} sequences")


if __name__ == "__main__":
    main()
