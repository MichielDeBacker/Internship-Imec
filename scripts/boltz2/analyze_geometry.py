from pathlib import Path
import csv
import json

import gemmi
import numpy as np

ROOT = Path.home() / "my_project"

PRED_ROOT = (
    ROOT
    / "results/boltz2_calibration/all"
    / "boltz_results_calibration_inputs/predictions"
)

CONF_SUMMARY = ROOT / "reports/boltz2/calibration_summary.tsv"
OUT = ROOT / "reports/boltz2/geometry_summary.tsv"

BACKBONE = {"N", "CA", "C", "O"}
CONTACT_CUTOFF = 5.0


def xyz(atom):
    return np.array(
        [atom.pos.x, atom.pos.y, atom.pos.z],
        dtype=float,
    )


def protein_residues(chain):
    out = []

    for res in chain:
        # Boltz protein residues should contain CA.
        names = {a.name.strip() for a in res}

        if "CA" in names:
            out.append(res)

    return out


def atoms_from_residue(res, mode="heavy"):
    coords = []

    for atom in res:
        name = atom.name.strip()

        if atom.element.name == "H":
            continue

        if mode == "backbone" and name not in BACKBONE:
            continue

        coords.append(xyz(atom))

    return coords


def chain_atoms(chain, mode="heavy"):
    coords = []

    for res in protein_residues(chain):
        coords.extend(atoms_from_residue(res, mode))

    if not coords:
        return np.empty((0, 3), dtype=float)

    return np.vstack(coords)


def minimum_distance(a, b):
    if len(a) == 0 or len(b) == 0:
        return float("nan")

    best = float("inf")

    # Chunk to avoid unnecessarily large temporary arrays.
    for start in range(0, len(a), 256):
        block = a[start:start + 256]
        delta = block[:, None, :] - b[None, :, :]
        d2 = np.sum(delta * delta, axis=2)
        best = min(best, float(np.sqrt(np.min(d2))))

    return best


def target_residue_contacts(target, binder, target_mode):
    binder_heavy = chain_atoms(binder, "heavy")
    count = 0

    for res in protein_residues(target):
        target_coords = atoms_from_residue(res, target_mode)

        if not target_coords:
            continue

        target_coords = np.vstack(target_coords)

        if minimum_distance(target_coords, binder_heavy) <= CONTACT_CUTOFF:
            count += 1

    return count


def read_confidences():
    data = {}

    with CONF_SUMMARY.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            data[row["name"]] = row

    return data


confidence = read_confidences()
rows = []

for cif in sorted(PRED_ROOT.glob("*/*_model_0.cif")):
    name = cif.parent.name

    structure = gemmi.read_structure(str(cif))
    model = structure[0]

    chains = {chain.name: chain for chain in model}

    missing = {"A", "B", "C"} - set(chains)

    if missing:
        raise RuntimeError(
            f"{name}: missing expected chains {sorted(missing)}; "
            f"found {sorted(chains)}"
        )

    A = chains["A"]
    B = chains["B"]
    C = chains["C"]

    A_bb = chain_atoms(A, "backbone")
    B_bb = chain_atoms(B, "backbone")
    C_bb = chain_atoms(C, "backbone")

    A_heavy = chain_atoms(A, "heavy")
    B_heavy = chain_atoms(B, "heavy")
    C_heavy = chain_atoms(C, "heavy")

    conf = confidence[name]

    ac = float(conf["AC_pair_iptm"])
    bc = float(conf["BC_pair_iptm"])

    rows.append({
        "name": name,
        "role": conf["role"],
        "len_A": len(protein_residues(A)),
        "len_B": len(protein_residues(B)),
        "len_C": len(protein_residues(C)),
        "AB_iptm": float(conf["AB_pair_iptm"]),
        "AC_iptm": ac,
        "BC_iptm": bc,
        "binder_iptm_mean": (ac + bc) / 2,
        "binder_iptm_min": min(ac, bc),
        "binder_iptm_delta": abs(ac - bc),
        "A_heavy_contacts": target_residue_contacts(A, C, "heavy"),
        "B_heavy_contacts": target_residue_contacts(B, C, "heavy"),
        "A_bbheavy_contacts": target_residue_contacts(A, C, "backbone"),
        "B_bbheavy_contacts": target_residue_contacts(B, C, "backbone"),
        "AC_bb_min": minimum_distance(A_bb, C_bb),
        "BC_bb_min": minimum_distance(B_bb, C_bb),
        "AC_heavy_min": minimum_distance(A_heavy, C_heavy),
        "BC_heavy_min": minimum_distance(B_heavy, C_heavy),
    })


rows.sort(
    key=lambda r: (
        -r["binder_iptm_min"],
        r["binder_iptm_delta"],
    )
)

fields = [
    "name",
    "role",
    "len_A",
    "len_B",
    "len_C",
    "AB_iptm",
    "AC_iptm",
    "BC_iptm",
    "binder_iptm_mean",
    "binder_iptm_min",
    "binder_iptm_delta",
    "A_heavy_contacts",
    "B_heavy_contacts",
    "A_bbheavy_contacts",
    "B_bbheavy_contacts",
    "AC_bb_min",
    "BC_bb_min",
    "AC_heavy_min",
    "BC_heavy_min",
]

with OUT.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

print(OUT.read_text())
