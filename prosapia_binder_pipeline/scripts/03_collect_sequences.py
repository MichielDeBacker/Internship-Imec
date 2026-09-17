#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from prosapia_binder_pipeline.common import (  # noqa: E402
    fasta_records,
    load_config,
    parse_kv_header,
    read_tsv,
    work_path,
    write_tsv,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    backbones = read_tsv(work_path(cfg, "01_backbones", "backbones.tsv"))
    nseq = int(cfg["sequencing"]["num_sequences_per_backbone"])

    out_rows = []
    missing = []

    for backbone in backbones:
        design_id = backbone["design_id"]
        fasta = work_path(cfg, "02_sequence", "runs", design_id, "seqs", f"{design_id}.fa")
        if not fasta.is_file():
            missing.append(design_id)
            continue

        records = fasta_records(fasta)
        if len(records) < 2:
            raise RuntimeError(f"No sampled ProteinMPNN sequences in {fasta}")

        sampled = records[1:]
        if len(sampled) != nseq:
            raise RuntimeError(
                f"{design_id}: expected {nseq} sampled sequences, found {len(sampled)}"
            )

        for ordinal, (header, sequence) in enumerate(sampled, start=1):
            meta = parse_kv_header(header)
            sample = int(meta.get("sample", ordinal))
            if "/" in sequence:
                raise RuntimeError(
                    f"{design_id}: expected only designed chain C in FASTA, got {sequence}"
                )
            if len(sequence) != int(backbone["chain_C_length"]):
                raise RuntimeError(
                    f"{design_id}: sequence length {len(sequence)} != chain_C_length {backbone['chain_C_length']}"
                )

            out_rows.append(
                {
                    "sequence_index": len(out_rows),
                    "sequence_id": f"{design_id}__s{sample:03d}",
                    "backbone_index": backbone["backbone_index"],
                    "design_id": design_id,
                    "role": backbone["role"],
                    "sample": sample,
                    "binder_sequence": sequence,
                    "binder_length": len(sequence),
                    "mpnn_temperature": meta.get("T", ""),
                    "mpnn_score": meta.get("score", ""),
                    "mpnn_global_score": meta.get("global_score", ""),
                    "mpnn_seq_recovery": meta.get("seq_recovery", ""),
                    "canonical_pdb": backbone["canonical_pdb"],
                    "target_template_pdb": backbone["target_template_pdb"],
                    "chain_A_sequence": backbone["chain_A_sequence"],
                    "chain_B_sequence": backbone["chain_B_sequence"],
                }
            )

    if missing:
        preview = "\n".join(missing[:20])
        raise RuntimeError(
            f"Missing sequence outputs for {len(missing)} backbones. First missing:\n{preview}"
        )

    expected = len(backbones) * nseq
    if len(out_rows) != expected:
        raise RuntimeError(f"Expected {expected} sequence rows, found {len(out_rows)}")

    fields = list(out_rows[0].keys())
    out = work_path(cfg, "02_sequence", "sequences.tsv")
    write_tsv(out, out_rows, fields)

    print(f"Backbones: {len(backbones)}")
    print(f"Sequences per backbone: {nseq}")
    print(f"Total binder sequences: {len(out_rows)}")
    print(f"Manifest: {out}")


if __name__ == "__main__":
    main()
