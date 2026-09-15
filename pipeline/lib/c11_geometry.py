#!/usr/bin/env python3

import argparse
import csv
import gzip
import io
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import biotite.structure.io.pdb as pdbio
import biotite.structure.io.pdbx as pdbx


BACKBONE = {"N", "CA", "C", "O"}


def load_structure(path):
    path = Path(path)
    low = path.name.lower()

    if low.endswith(".cif.gz") or low.endswith(".mmcif.gz"):
        with gzip.open(path, "rt") as fh:
            f = pdbx.CIFFile.read(io.StringIO(fh.read()))
        return pdbx.get_structure(f, model=1)

    if low.endswith(".cif") or low.endswith(".mmcif"):
        f = pdbx.CIFFile.read(path)
        return pdbx.get_structure(f, model=1)

    if low.endswith(".pdb"):
        f = pdbio.PDBFile.read(path)
        return f.get_structure(model=1)

    raise ValueError(f"Unsupported structure: {path}")


def chain_ids(a):
    out = []
    for c in a.chain_id:
        c = str(c)
        if c not in out:
            out.append(c)
    return out


def chain(a, cid):
    return a[np.asarray(a.chain_id, dtype=str) == str(cid)]


def ordered_ca(a):
    """
    One CA coordinate per residue, preserving structure order.
    """
    xyz = []
    resid = []
    seen = set()

    for i in range(len(a)):
        if str(a.atom_name[i]) != "CA":
            continue

        key = (
            str(a.chain_id[i]),
            int(a.res_id[i]),
        )

        if key in seen:
            continue

        seen.add(key)
        xyz.append(np.asarray(a.coord[i], dtype=float))
        resid.append(int(a.res_id[i]))

    return np.asarray(xyz, dtype=float), resid


def backbone_coords(a):
    mask = np.asarray(
        [str(x) in BACKBONE for x in a.atom_name],
        dtype=bool,
    )
    return np.asarray(a.coord[mask], dtype=float)


def heavy_coords(a):
    # Hydrogen is normally absent anyway, but make it explicit.
    if hasattr(a, "element"):
        mask = np.asarray(
            [str(x).upper() != "H" for x in a.element],
            dtype=bool,
        )
        return np.asarray(a.coord[mask], dtype=float), np.where(mask)[0]

    idx = np.arange(len(a))
    return np.asarray(a.coord, dtype=float), idx


def kabsch(mobile, target):
    mobile = np.asarray(mobile, dtype=float)
    target = np.asarray(target, dtype=float)

    if mobile.shape != target.shape:
        raise ValueError(
            f"Kabsch shape mismatch: {mobile.shape} vs {target.shape}"
        )

    if len(mobile) < 3:
        raise ValueError("Need at least 3 points for alignment")

    cm = mobile.mean(axis=0)
    ct = target.mean(axis=0)

    x = mobile - cm
    y = target - ct

    h = x.T @ y
    u, s, vt = np.linalg.svd(h)

    r = vt.T @ u.T

    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T

    # Coordinates are row vectors.
    t = ct - cm @ r.T

    fitted = mobile @ r.T + t

    rmsd = float(
        np.sqrt(
            np.mean(
                np.sum(
                    (fitted - target) ** 2,
                    axis=1,
                )
            )
        )
    )

    return r, t, rmsd


def transform_coords(coords, r, t):
    return np.asarray(coords, dtype=float) @ r.T + t


def transformed_array(a, r, t):
    b = a.copy()
    b.coord = transform_coords(a.coord, r, t)
    return b


def min_distance(a, b):
    if len(a) == 0 or len(b) == 0:
        return math.inf

    tree = cKDTree(b)
    d, _ = tree.query(a, k=1)

    return float(np.min(d))


def contacted_binder_residues(binder, target, cutoff=5.0):
    bxyz, bind_idx = heavy_coords(binder)
    txyz, _ = heavy_coords(target)

    if len(bxyz) == 0 or len(txyz) == 0:
        return set()

    tree = cKDTree(txyz)

    hits = tree.query_ball_point(
        bxyz,
        r=float(cutoff),
    )

    residues = set()

    for local_i, nbrs in enumerate(hits):
        if not nbrs:
            continue

        original_i = int(bind_idx[local_i])

        residues.add(
            (
                str(binder.chain_id[original_i]),
                int(binder.res_id[original_i]),
            )
        )

    return residues


def detect_target_binder(design, native_chain):
    """
    Parent RFD3 structures contain target + binder.
    Identify target as the chain whose CA count most closely matches
    a native SlyB chain; the other chain is the binder.
    """
    ids = chain_ids(design)

    if len(ids) != 2:
        raise ValueError(
            f"Expected 2 chains in parent pose, found {ids}"
        )

    native_ca, _ = ordered_ca(native_chain)
    n_native = len(native_ca)

    stats = []

    for cid in ids:
        ca, _ = ordered_ca(chain(design, cid))
        stats.append(
            (
                abs(len(ca) - n_native),
                -len(ca),
                cid,
                len(ca),
            )
        )

    stats.sort()

    target_id = stats[0][2]
    binder_id = [x for x in ids if x != target_id][0]

    target_ca, _ = ordered_ca(chain(design, target_id))

    if len(target_ca) != n_native:
        raise ValueError(
            "Target/native CA count mismatch: "
            f"design chain {target_id}={len(target_ca)}, "
            f"native={n_native}"
        )

    return target_id, binder_id


def identify_partner_for_transform(
    ring,
    ring_ids,
    ref_b_ca,
    r,
    t,
    primary_id,
):
    """
    Apply A->current transform to native reference B and find which
    real C11 chain it lands on.
    """
    predicted = transform_coords(ref_b_ca, r, t)

    candidates = []

    for cid in ring_ids:
        if cid == primary_id:
            continue

        ca, _ = ordered_ca(chain(ring, cid))

        if ca.shape != predicted.shape:
            continue

        rmsd = float(
            np.sqrt(
                np.mean(
                    np.sum(
                        (predicted - ca) ** 2,
                        axis=1,
                    )
                )
            )
        )

        candidates.append((rmsd, cid))

    if not candidates:
        raise ValueError(
            f"No valid partner for ring chain {primary_id}"
        )

    candidates.sort()

    return candidates[0][1], candidates[0][0]


def evaluate(
    ring_path,
    parent_path,
    monomer_path=None,
    name=None,
):
    ring = load_structure(ring_path)
    parent = load_structure(parent_path)

    ring_ids = chain_ids(ring)

    if len(ring_ids) != 11:
        raise ValueError(
            f"Expected C11 ring with 11 chains, found {ring_ids}"
        )

    # Project convention/native reference pair.
    if "A" not in ring_ids or "B" not in ring_ids:
        raise ValueError(
            f"Native ring must contain chains A and B: {ring_ids}"
        )

    ref_a = chain(ring, "A")
    ref_b = chain(ring, "B")

    ref_a_ca, _ = ordered_ca(ref_a)
    ref_b_ca, _ = ordered_ca(ref_b)

    target_id, binder_id = detect_target_binder(
        parent,
        ref_a,
    )

    parent_target = chain(parent, target_id)
    parent_binder = chain(parent, binder_id)

    parent_target_ca, _ = ordered_ca(parent_target)
    parent_binder_ca, _ = ordered_ca(parent_binder)

    # Parent target -> native reference chain A.
    r_parent, t_parent, target_fit_rmsd = kabsch(
        parent_target_ca,
        ref_a_ca,
    )

    binder_source = parent_binder
    monomer_fit_rmsd = None

    if monomer_path:
        monomer = load_structure(monomer_path)
        mids = chain_ids(monomer)

        if len(mids) != 1:
            raise ValueError(
                f"RF3 monomer expected 1 chain, found {mids}"
            )

        mono = chain(monomer, mids[0])
        mono_ca, _ = ordered_ca(mono)

        if mono_ca.shape != parent_binder_ca.shape:
            raise ValueError(
                "RF3 monomer/backbone CA mismatch: "
                f"{len(mono_ca)} vs {len(parent_binder_ca)}"
            )

        # RF3 monomer -> parent binder pose.
        r_mono, t_mono, monomer_fit_rmsd = kabsch(
            mono_ca,
            parent_binder_ca,
        )

        binder_source = transformed_array(
            mono,
            r_mono,
            t_mono,
        )

    # Parent design frame -> native ring-A frame.
    binder_ref = transformed_array(
        binder_source,
        r_parent,
        t_parent,
    )

    placements = []
    partner_rows = []

    for cid in ring_ids:
        current = chain(ring, cid)
        current_ca, _ = ordered_ca(current)

        if current_ca.shape != ref_a_ca.shape:
            raise ValueError(
                f"Native chain {cid} does not match chain A"
            )

        r_sym, t_sym, a_residual = kabsch(
            ref_a_ca,
            current_ca,
        )

        partner_id, neighbor_residual = (
            identify_partner_for_transform(
                ring,
                ring_ids,
                ref_b_ca,
                r_sym,
                t_sym,
                cid,
            )
        )

        placed_binder = transformed_array(
            binder_ref,
            r_sym,
            t_sym,
        )

        target1 = current
        target2 = chain(ring, partner_id)

        contacts1 = contacted_binder_residues(
            placed_binder,
            target1,
            cutoff=5.0,
        )

        contacts2 = contacted_binder_residues(
            placed_binder,
            target2,
            cutoff=5.0,
        )

        target_min = min_distance(
            backbone_coords(placed_binder),
            backbone_coords(ring),
        )

        placements.append(
            {
                "primary_chain": cid,
                "partner_chain": partner_id,
                "binder": placed_binder,
                "target_min_bb_A": target_min,
                "contact_res_partner1": len(contacts1),
                "contact_res_partner2": len(contacts2),
                "A_fit_rmsd_A": a_residual,
                "neighbor_residual_A": neighbor_residual,
            }
        )

        partner_rows.append(
            (
                cid,
                partner_id,
                neighbor_residual,
            )
        )

    # All 55 binder-copy pairs.
    binder_pair_min = math.inf

    for i in range(len(placements)):
        bi = backbone_coords(
            placements[i]["binder"]
        )

        for j in range(i + 1, len(placements)):
            bj = backbone_coords(
                placements[j]["binder"]
            )

            d = min_distance(bi, bj)

            binder_pair_min = min(
                binder_pair_min,
                d,
            )

    min_target_bb = min(
        x["target_min_bb_A"]
        for x in placements
    )

    min_contacts_1 = min(
        x["contact_res_partner1"]
        for x in placements
    )

    min_contacts_2 = min(
        x["contact_res_partner2"]
        for x in placements
    )

    max_neighbor_residual = max(
        x["neighbor_residual_A"]
        for x in placements
    )

    # Reuse-first project thresholds from current SlyB plan.
    thresholds = {
        "target_fit_rmsd_A_max": 2.0,
        "neighbor_residual_A_max": 1.0,
        "hard_bb_distance_A_min": 1.8,
        "distinct_contact_residues_per_partner_min": 3,
    }

    pass_geometry = bool(
        target_fit_rmsd <= thresholds["target_fit_rmsd_A_max"]
        and max_neighbor_residual
        <= thresholds["neighbor_residual_A_max"]
        and min_target_bb
        >= thresholds["hard_bb_distance_A_min"]
        and binder_pair_min
        >= thresholds["hard_bb_distance_A_min"]
        and min_contacts_1
        >= thresholds[
            "distinct_contact_residues_per_partner_min"
        ]
        and min_contacts_2
        >= thresholds[
            "distinct_contact_residues_per_partner_min"
        ]
    )

    result = {
        "name": name or Path(parent_path).stem,
        "parent": str(parent_path),
        "monomer": (
            str(monomer_path)
            if monomer_path
            else None
        ),
        "target_chain_in_parent": target_id,
        "binder_chain_in_parent": binder_id,
        "target_CA_count": len(parent_target_ca),
        "binder_CA_count": len(parent_binder_ca),
        "target_fit_CA_RMSD_A": target_fit_rmsd,
        "monomer_to_parent_CA_RMSD_A": monomer_fit_rmsd,
        "max_neighbor_placement_residual_A": (
            max_neighbor_residual
        ),
        "min_binder_to_ring_backbone_A": (
            min_target_bb
        ),
        "min_binder_to_binder_backbone_A": (
            binder_pair_min
        ),
        "min_contacted_binder_res_partner1": (
            min_contacts_1
        ),
        "min_contacted_binder_res_partner2": (
            min_contacts_2
        ),
        "placements": [
            {
                k: v
                for k, v in x.items()
                if k != "binder"
            }
            for x in placements
        ],
        "thresholds": thresholds,
        "PASS": pass_geometry,
    }

    return result


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--ring",
        required=True,
    )

    ap.add_argument(
        "--parent",
        required=True,
    )

    ap.add_argument(
        "--monomer",
    )

    ap.add_argument(
        "--name",
    )

    ap.add_argument(
        "--json-out",
    )

    args = ap.parse_args()

    result = evaluate(
        args.ring,
        args.parent,
        args.monomer,
        args.name,
    )

    text = json.dumps(
        result,
        indent=2,
    )

    print(text)

    if args.json_out:
        p = Path(args.json_out)
        p.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        p.write_text(text + "\n")


if __name__ == "__main__":
    main()
