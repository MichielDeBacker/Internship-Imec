#!/usr/bin/env python3

from pathlib import Path
from collections import Counter
import csv
import gzip
import io
import math

import numpy as np
from scipy.spatial import cKDTree

import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx


BASE = Path("/scratch/brussel/vo/000/bvo00014/vsc39230")
REPO = BASE / "repos/Internship-Imec"

TARGET = REPO / "work/vertical96_150/target_AB_historical96.pdb"
RING = BASE / "slyb/inputs/target_ring_C11.pdb"
OUTROOT = BASE / "slyb/out"

CSV_OUT = REPO / "work/vertical96_150/vertical96_screen.csv"
SURVIVORS_OUT = REPO / "work/vertical96_150/vertical96_survivors.txt"

TARGET_RMSD_MAX = 2.0
NEIGHBOR_RESIDUAL_MAX = 1.0
CONTACT_CUTOFF = 5.0
MIN_CONTACT_RESIDUES = 3
VERTICAL_ANGLE_MAX = 20.0
MIN_BB_DISTANCE = 1.8

BB_NAMES = {"N", "CA", "C", "O"}


def load_structure(path):
    path = Path(path)

    if path.name.endswith(".cif.gz"):
        with gzip.open(path, "rt") as f:
            cif = pdbx.CIFFile.read(io.StringIO(f.read()))
        return pdbx.get_structure(cif, model=1)

    if path.suffix == ".cif":
        cif = pdbx.CIFFile.read(str(path))
        return pdbx.get_structure(cif, model=1)

    if path.suffix == ".pdb":
        pf = pdb.PDBFile.read(str(path))
        return pf.get_structure(model=1)

    raise ValueError(f"unsupported file: {path}")


def chain_ids(arr):
    return list(dict.fromkeys(map(str, arr.chain_id)))


def chain(arr, cid):
    return arr[np.asarray(arr.chain_id, dtype=str) == cid]


def ca(arr):
    return arr[np.asarray(arr.atom_name, dtype=str) == "CA"]


def ca_xyz(arr):
    return np.asarray(ca(arr).coord, dtype=float)


def n_ca(arr):
    return len(ca(arr))


def historical96_ca(arr):
    rid = np.asarray(arr.res_id)
    name = np.asarray(arr.atom_name, dtype=str)

    mask = (
        (name == "CA")
        & (
            ((rid >= 18) & (rid <= 60))
            | ((rid >= 103) & (rid <= 155))
        )
    )

    out = np.asarray(arr.coord[mask], dtype=float)

    if len(out) != 96:
        raise ValueError(
            f"expected 96 historical96 CA atoms, got {len(out)}"
        )

    return out


def kabsch(P, Q):
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)

    if P.shape != Q.shape:
        raise ValueError(f"Kabsch shape mismatch: {P.shape} vs {Q.shape}")

    pc = P.mean(axis=0)
    qc = Q.mean(axis=0)

    X = P - pc
    Y = Q - qc

    U, s, Vt = np.linalg.svd(X.T @ Y)
    R = Vt.T @ U.T

    if np.linalg.det(R) < 0:
        Vt[-1] *= -1
        R = Vt.T @ U.T

    fitted = X @ R.T + qc
    rmsd = np.sqrt(np.mean(np.sum((fitted - Q) ** 2, axis=1)))

    return R, pc, qc, float(rmsd)


def transform(xyz, fit):
    R, pc, qc, _ = fit
    xyz = np.asarray(xyz, dtype=float)
    return (xyz - pc) @ R.T + qc


def direct_rmsd(P, Q):
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)

    if P.shape != Q.shape:
        return float("inf")

    return float(
        np.sqrt(np.mean(np.sum((P - Q) ** 2, axis=1)))
    )


def ring_geometry(ring):
    cids = chain_ids(ring)

    centroids = np.array([
        ca_xyz(chain(ring, c)).mean(axis=0)
        for c in cids
    ])

    center = centroids.mean(axis=0)

    X = centroids - center
    vals, vecs = np.linalg.eigh(X.T @ X)

    # Ring normal = direction with least centroid spread.
    axis = vecs[:, np.argmin(vals)]
    axis = axis / np.linalg.norm(axis)

    e1 = centroids[0] - center
    e1 = e1 - axis * np.dot(e1, axis)
    e1 = e1 / np.linalg.norm(e1)

    e2 = np.cross(axis, e1)

    angles = [
        math.atan2(
            np.dot(c - center, e2),
            np.dot(c - center, e1),
        )
        for c in centroids
    ]

    order_idx = np.argsort(angles)
    ordered = [cids[i] for i in order_idx]

    return center, axis, ordered


def map_reference_pair_to_ring(reference, ring, ring_order):
    rc = chain_ids(reference)

    if len(rc) != 2:
        raise ValueError(
            f"reference target should have 2 chains, found {rc}"
        )

    P1 = ca_xyz(chain(reference, rc[0]))
    P2 = ca_xyz(chain(reference, rc[1]))

    if len(P1) != 96 or len(P2) != 96:
        raise ValueError(
            f"reference chain lengths are {len(P1)}, {len(P2)}, expected 96/96"
        )

    P = np.concatenate([P1, P2], axis=0)

    best = None
    n = len(ring_order)

    for i in range(n):
        a = ring_order[i]
        b = ring_order[(i + 1) % n]

        for pair in [(a, b), (b, a)]:
            Q1 = historical96_ca(chain(ring, pair[0]))
            Q2 = historical96_ca(chain(ring, pair[1]))
            Q = np.concatenate([Q1, Q2], axis=0)

            fit = kabsch(P, Q)

            candidate = {
                "rmsd": fit[3],
                "fit": fit,
                "reference_chains": tuple(rc),
                "ring_pair": pair,
            }

            if best is None or candidate["rmsd"] < best["rmsd"]:
                best = candidate

    return best


def map_candidate_to_reference(arr, reference):
    cids = chain_ids(arr)

    if len(cids) != 3:
        raise ValueError(
            f"expected 3 chains in RFD3 output, found {cids}"
        )

    # Shortest chain is binder; two long chains are target.
    ordered = sorted(cids, key=lambda c: n_ca(chain(arr, c)))

    binder_cid = ordered[0]
    target_cids = ordered[1:]

    rc = chain_ids(reference)

    ref1 = ca_xyz(chain(reference, rc[0]))
    ref2 = ca_xyz(chain(reference, rc[1]))

    best = None

    for cand_order in [
        target_cids,
        list(reversed(target_cids)),
    ]:
        c1 = ca_xyz(chain(arr, cand_order[0]))
        c2 = ca_xyz(chain(arr, cand_order[1]))

        if c1.shape != ref1.shape or c2.shape != ref2.shape:
            continue

        P = np.concatenate([c1, c2], axis=0)
        Q = np.concatenate([ref1, ref2], axis=0)

        pair_fit = kabsch(P, Q)

        # Anchor target 1, then inspect placement of target 2.
        fit1 = kabsch(c1, ref1)
        residual_1_to_2 = direct_rmsd(
            transform(c2, fit1),
            ref2,
        )

        # Reverse anchor as a symmetric check.
        fit2 = kabsch(c2, ref2)
        residual_2_to_1 = direct_rmsd(
            transform(c1, fit2),
            ref1,
        )

        neighbor_residual = max(
            residual_1_to_2,
            residual_2_to_1,
        )

        candidate = {
            "binder_cid": binder_cid,
            "target_cids": tuple(cand_order),
            "pair_fit": pair_fit,
            "target_rmsd": pair_fit[3],
            "neighbor_residual": neighbor_residual,
            "residual_1_to_2": residual_1_to_2,
            "residual_2_to_1": residual_2_to_1,
        }

        if best is None or candidate["target_rmsd"] < best["target_rmsd"]:
            best = candidate

    if best is None:
        raise ValueError("could not map generated target pair to reference")

    return best


def heavy_mask(arr):
    elem = np.char.upper(np.asarray(arr.element, dtype=str))
    names = np.char.upper(np.asarray(arr.atom_name, dtype=str))

    return (elem != "H") & (~np.char.startswith(names, "H"))


def backbone_mask(arr):
    names = np.asarray(arr.atom_name, dtype=str)
    return np.isin(names, list(BB_NAMES))


def contact_residue_count(
    binder_arr,
    binder_xyz,
    target_arr,
    cutoff=CONTACT_CUTOFF,
):
    bmask = heavy_mask(binder_arr)
    tmask = heavy_mask(target_arr)

    bx = np.asarray(binder_xyz, dtype=float)[bmask]
    tx = np.asarray(target_arr.coord, dtype=float)[tmask]

    if len(bx) == 0 or len(tx) == 0:
        return 0

    tree = cKDTree(tx)
    distances, _ = tree.query(bx, k=1)

    close = distances <= cutoff

    res_ids = np.asarray(binder_arr.res_id)[bmask]

    return len(set(map(int, res_ids[close])))


def principal_axis(xyz):
    xyz = np.asarray(xyz, dtype=float)
    X = xyz - xyz.mean(axis=0)

    vals, vecs = np.linalg.eigh(X.T @ X)

    axis = vecs[:, np.argmax(vals)]
    return axis / np.linalg.norm(axis)


def angle_deg(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)

    # Axis has no direction, so +v and -v are equivalent.
    c = np.clip(abs(np.dot(a, b)), 0.0, 1.0)

    return float(np.degrees(np.arccos(c)))


def rotation_matrix(axis, theta):
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)

    x, y, z = axis
    c = math.cos(theta)
    s = math.sin(theta)
    C = 1.0 - c

    return np.array([
        [c + x*x*C,     x*y*C - z*s, x*z*C + y*s],
        [y*x*C + z*s,   c + y*y*C,   y*z*C - x*s],
        [z*x*C - y*s,   z*y*C + x*s, c + z*z*C],
    ])


def rotate_about_axis(xyz, center, axis, theta):
    R = rotation_matrix(axis, theta)
    xyz = np.asarray(xyz, dtype=float)

    return (xyz - center) @ R.T + center


def minimum_copy_backbone_distance(copy_bb_xyz):
    best = float("inf")

    for i in range(len(copy_bb_xyz)):
        tree = cKDTree(copy_bb_xyz[i])

        for j in range(i + 1, len(copy_bb_xyz)):
            d, _ = tree.query(copy_bb_xyz[j], k=1)
            best = min(best, float(np.min(d)))

    return best


def screen_one(
    path,
    reference,
    ring,
    ring_center,
    ring_axis,
    ring_order,
    ref_to_ring,
):
    arr = load_structure(path)

    mapping = map_candidate_to_reference(arr, reference)

    binder = chain(arr, mapping["binder_cid"])

    # Candidate -> canonical historical96 pair.
    binder_ref = transform(
        np.asarray(binder.coord, dtype=float),
        mapping["pair_fit"],
    )

    binder_ref_ca = transform(
        ca_xyz(binder),
        mapping["pair_fit"],
    )

    # Canonical historical96 pair -> actual native C11 pair.
    binder_native = transform(
        binder_ref,
        ref_to_ring["fit"],
    )

    binder_native_ca = transform(
        binder_ref_ca,
        ref_to_ring["fit"],
    )

    vertical_angle = angle_deg(
        principal_axis(binder_native_ca),
        ring_axis,
    )

    pair0 = ref_to_ring["ring_pair"]

    ia = ring_order.index(pair0[0])
    ib = ring_order.index(pair0[1])

    n = len(ring_order)

    bbmask = backbone_mask(binder)

    copy_bb = []
    copy_contacts = []

    for k in range(n):
        theta = 2.0 * math.pi * k / n

        xyz = rotate_about_axis(
            binder_native,
            ring_center,
            ring_axis,
            theta,
        )

        copy_bb.append(xyz[bbmask])

        partner1 = ring_order[(ia + k) % n]
        partner2 = ring_order[(ib + k) % n]

        c1 = contact_residue_count(
            binder,
            xyz,
            chain(ring, partner1),
        )

        c2 = contact_residue_count(
            binder,
            xyz,
            chain(ring, partner2),
        )

        copy_contacts.append(
            (partner1, partner2, c1, c2)
        )

    min_contacts_1 = min(x[2] for x in copy_contacts)
    min_contacts_2 = min(x[3] for x in copy_contacts)

    contacts_ok = all(
        x[2] >= MIN_CONTACT_RESIDUES
        and x[3] >= MIN_CONTACT_RESIDUES
        for x in copy_contacts
    )

    min_bb_distance = minimum_copy_backbone_distance(copy_bb)

    target_ok = (
        mapping["target_rmsd"]
        <= TARGET_RMSD_MAX
    )

    neighbor_ok = (
        mapping["neighbor_residual"]
        <= NEIGHBOR_RESIDUAL_MAX
    )

    vertical_ok = (
        vertical_angle
        <= VERTICAL_ANGLE_MAX
    )

    clash_ok = (
        min_bb_distance
        >= MIN_BB_DISTANCE
    )

    reasons = []

    if not target_ok:
        reasons.append("target_rmsd")

    if not neighbor_ok:
        reasons.append("neighbor_placement")

    if not contacts_ok:
        reasons.append("C11_contacts")

    if not vertical_ok:
        reasons.append("vertical_angle")

    if not clash_ok:
        reasons.append("binder_binder_clash")

    passed = (
        target_ok
        and neighbor_ok
        and contacts_ok
        and vertical_ok
        and clash_ok
    )

    return {
        "file": str(path),
        "binder_chain": mapping["binder_cid"],
        "binder_length": n_ca(binder),
        "target_chain_1": mapping["target_cids"][0],
        "target_chain_2": mapping["target_cids"][1],

        "target_CA_RMSD_A": round(
            mapping["target_rmsd"], 3
        ),

        "neighbor_residual_A": round(
            mapping["neighbor_residual"], 3
        ),

        "neighbor_1_to_2_A": round(
            mapping["residual_1_to_2"], 3
        ),

        "neighbor_2_to_1_A": round(
            mapping["residual_2_to_1"], 3
        ),

        "vertical_angle_deg": round(
            vertical_angle, 2
        ),

        "C11_min_contacts_partner1": min_contacts_1,
        "C11_min_contacts_partner2": min_contacts_2,
        "C11_all_copies_contact_ok": contacts_ok,

        "min_binder_binder_BB_distance_A": round(
            min_bb_distance, 3
        ),

        "status": (
            "NEEDS_VISUAL_REVIEW"
            if passed
            else "REJECT_OR_REPAIR"
        ),

        "reasons": ";".join(reasons),
    }


def find_backbones():
    files = []

    for jobdir in OUTROOT.glob("rfd3_vertical96_150_*"):
        files.extend(
            jobdir.rglob("*_model_*.cif.gz")
        )

    # Remove accidental duplicate path entries.
    files = sorted(set(p.resolve() for p in files))

    return files


def main():
    reference = load_structure(TARGET)
    ring = load_structure(RING)

    ring_center, ring_axis, ring_order = ring_geometry(ring)

    if len(ring_order) != 11:
        raise RuntimeError(
            f"expected C11 ring, found {len(ring_order)} chains"
        )

    ref_to_ring = map_reference_pair_to_ring(
        reference,
        ring,
        ring_order,
    )

    print("=== CANONICAL TARGET -> NATIVE C11 ===")
    print("ring pair:", ref_to_ring["ring_pair"])
    print(
        "mapping RMSD:",
        f"{ref_to_ring['rmsd']:.3f} A",
    )

    if ref_to_ring["rmsd"] > 2.0:
        raise RuntimeError(
            "canonical target does not map cleanly to native C11; "
            "stop before screening"
        )

    files = find_backbones()

    print()
    print("=== INPUT ===")
    print("backbones found:", len(files))

    if len(files) != 150:
        print(
            f"WARNING: expected 150 files, found {len(files)}"
        )

    rows = []

    for i, path in enumerate(files, 1):
        try:
            row = screen_one(
                path,
                reference,
                ring,
                ring_center,
                ring_axis,
                ring_order,
                ref_to_ring,
            )

        except Exception as e:
            row = {
                "file": str(path),
                "binder_chain": "",
                "binder_length": "",
                "target_chain_1": "",
                "target_chain_2": "",
                "target_CA_RMSD_A": "",
                "neighbor_residual_A": "",
                "neighbor_1_to_2_A": "",
                "neighbor_2_to_1_A": "",
                "vertical_angle_deg": "",
                "C11_min_contacts_partner1": "",
                "C11_min_contacts_partner2": "",
                "C11_all_copies_contact_ok": False,
                "min_binder_binder_BB_distance_A": "",
                "status": "REJECT_OR_REPAIR",
                "reasons": f"geometry_error:{type(e).__name__}:{e}",
            }

        rows.append(row)

        print(
            f"[{i:3d}/{len(files)}] "
            f"{path.name} -> "
            f"{row['status']} "
            f"{row['reasons']}"
        )

    if not rows:
        raise RuntimeError("no RFD3 backbone files found")

    CSV_OUT.parent.mkdir(parents=True, exist_ok=True)

    with open(CSV_OUT, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)

    survivors = [
        r for r in rows
        if r["status"] == "NEEDS_VISUAL_REVIEW"
    ]

    with open(SURVIVORS_OUT, "w") as f:
        for r in survivors:
            f.write(r["file"] + "\n")

    reject_reasons = Counter()

    for row in rows:
        if row["reasons"]:
            for reason in row["reasons"].split(";"):
                reject_reasons[reason] += 1

    print()
    print("========================================")
    print("SCREEN SUMMARY")
    print("========================================")
    print("screened:", len(rows))
    print("NEEDS_VISUAL_REVIEW:", len(survivors))
    print("REJECT_OR_REPAIR:", len(rows) - len(survivors))

    print()
    print("Rejection reasons:")
    for reason, count in reject_reasons.most_common():
        print(f"  {reason}: {count}")

    print()
    print("CSV:", CSV_OUT)
    print("survivor list:", SURVIVORS_OUT)


if __name__ == "__main__":
    main()
