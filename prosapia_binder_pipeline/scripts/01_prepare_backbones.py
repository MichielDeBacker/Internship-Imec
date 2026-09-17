#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
import glob
import gzip
import sys
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.SeqUtils import seq1

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from prosapia_binder_pipeline.common import load_config, safe_id, work_path, write_tsv  # noqa: E402


class TargetOnly(Select):
    def accept_chain(self, chain):
        return chain.id in {"A", "B"}


def ca_residues(chain):
    return [res for res in chain if res.resname != "HOH" and "CA" in res]


def chain_sequence(chain) -> str:
    letters = []
    for residue in ca_residues(chain):
        letters.append(seq1(residue.resname, undef_code="X"))
    return "".join(letters)


def classify(model, cfg):
    bcfg = cfg["backbones"]
    info = [(chain.id, len(ca_residues(chain)), chain) for chain in model]
    info = [x for x in info if x[1] > 0]

    binders = [
        x for x in info
        if bcfg["binder_length_min"] <= x[1] <= bcfg["binder_length_max"]
    ]
    targets = [
        x for x in info
        if bcfg["target_length_min"] <= x[1] <= bcfg["target_length_max"]
    ]

    if len(binders) != 1 or len(targets) != 2:
        observed = [(x[0], x[1]) for x in info]
        raise RuntimeError(
            "Could not uniquely classify one binder and two target chains. "
            f"Observed chain lengths: {observed}"
        )

    # The two target chains are symmetry-equivalent. Keep deterministic ordering.
    targets.sort(key=lambda x: str(x[0]))
    return targets[0], targets[1], binders[0]


def canonical_structure(name, target_a, target_b, binder):
    structure = Structure(name)
    model = Model(0)
    structure.add(model)
    for new_id, entry in zip(("A", "B", "C"), (target_a, target_b, binder)):
        chain = copy.deepcopy(entry[2])
        chain.id = new_id
        model.add(chain)
    return structure


def make_design_id(path: Path) -> str:
    filename = path.name
    if filename.endswith(".cif.gz"):
        filename = filename[:-7]
    elif filename.endswith(".cif"):
        filename = filename[:-4]
    return safe_id(f"{path.parent.name}__{filename}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    expected = int(cfg["project"]["expected_backbones"])
    paths = [Path(p) for p in sorted(glob.glob(cfg["paths"]["backbone_glob"]))]

    if len(paths) != expected:
        raise RuntimeError(
            f"Expected exactly {expected} RFD3 backbones, found {len(paths)}. "
            "Refusing to continue so a partial/duplicate set is not silently sequenced."
        )

    input_dir = work_path(cfg, "01_backbones", "canonical_pdb")
    template_dir = work_path(cfg, "01_backbones", "target_templates")
    input_dir.mkdir(parents=True, exist_ok=True)
    template_dir.mkdir(parents=True, exist_ok=True)

    parser = MMCIFParser(QUIET=True)
    rows = []
    seen = set()

    for index, source in enumerate(paths):
        design_id = make_design_id(source)
        if design_id in seen:
            raise RuntimeError(f"Duplicate design_id after normalization: {design_id}")
        seen.add(design_id)

        with gzip.open(source, "rt") as handle:
            structure = parser.get_structure(design_id, handle)

        models = list(structure.get_models())
        if len(models) != 1:
            raise RuntimeError(f"{source}: expected one model, found {len(models)}")

        target_a, target_b, binder = classify(models[0], cfg)
        canonical = canonical_structure(design_id, target_a, target_b, binder)

        canonical_pdb = input_dir / f"{design_id}.pdb"
        target_template = template_dir / f"{design_id}__AB_template.pdb"

        io = PDBIO()
        io.set_structure(canonical)
        io.save(str(canonical_pdb))
        io.save(str(target_template), TargetOnly())

        can_model = next(canonical.get_models())
        chain_a = can_model["A"]
        chain_b = can_model["B"]
        chain_c = can_model["C"]

        row = {
            "backbone_index": index,
            "design_id": design_id,
            "role": cfg["project"].get("role", "candidate"),
            "source_cif_gz": str(source),
            "canonical_pdb": str(canonical_pdb),
            "target_template_pdb": str(target_template),
            "original_target_a_chain": target_a[0],
            "original_target_b_chain": target_b[0],
            "original_binder_chain": binder[0],
            "chain_A_length": len(ca_residues(chain_a)),
            "chain_B_length": len(ca_residues(chain_b)),
            "chain_C_length": len(ca_residues(chain_c)),
            "chain_A_sequence": chain_sequence(chain_a),
            "chain_B_sequence": chain_sequence(chain_b),
        }
        rows.append(row)
        print(
            f"[{index + 1:03d}/{expected}] {design_id} "
            f"A={row['chain_A_length']} B={row['chain_B_length']} C={row['chain_C_length']}"
        )

    fields = list(rows[0].keys())
    manifest = work_path(cfg, "01_backbones", "backbones.tsv")
    write_tsv(manifest, rows, fields)
    print(f"\nWrote {len(rows)} canonical complexes and AB templates")
    print(f"Manifest: {manifest}")


if __name__ == "__main__":
    main()
