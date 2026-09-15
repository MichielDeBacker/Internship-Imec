#!/usr/bin/env python3

import argparse
import csv
import gzip
import io
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import biotite.structure as struc
import biotite.structure.io.pdb as pdbio
import biotite.structure.io.pdbx as pdbx


# Current implemented thresholds from the existing SlyB filter
CLASH_A = 2.2
MIN_NBR_A = 2.2
COMFORT_A = 3.4
CONTACT_A = 5.0
MIN_CONTACTS = 20

# Target geometry should be preserved after rigid-body alignment
MAX_TARGET_CA_RMSD_A = 1.5


def load_structure(path):
    path = str(path)
    low = path.lower()

    if low.endswith((".cif", ".cif.gz", ".mmcif", ".mmcif.gz")):
        if low.endswith(".gz"):
            with gzip.open(path, "rt") as fh:
                f = pdbx.CIFFile.read(io.StringIO(fh.read()))
        else:
            f = pdbx.CIFFile.read(path)

        a = pdbx.get_structure(f, model=1)
    else:
        a = pdbio.PDBFile.read(path).get_structure(model=1)

    return a[(a.element != "H") & struc.filter_amino_acids(a)]


def parse_label(x):
    m = re.fullmatch(r"(.+?)(-?\d+)", str(x))
    if not m:
        raise ValueError(f"cannot parse residue label: {x}")
    return m.group(1), int(m.group(2))


def sidecar_for(cif):
    s = str(cif)
    if not s.endswith(".cif.gz"):
        raise ValueError(f"expected .cif.gz: {cif}")
    return Path(s[:-7] + ".json")


def kabsch(moving, fixed):
    moving = np.asarray(moving, dtype=float)
    fixed = np.asarray(fixed, dtype=float)

    mc = moving.mean(axis=0)
    fc = fixed.mean(axis=0)

    P = moving - mc
    Q = fixed - fc

    U, S, Vt = np.linalg.svd(P.T @ Q)
    R = U @ Vt

    if np.linalg.det(R) < 0:
        U[:, -1] *= -1
        R = U @ Vt

    aligned = P @ R + fc
    rmsd = float(
        np.sqrt(np.mean(np.sum((aligned - fixed) ** 2, axis=1)))
    )

    return R, mc, fc, rmsd


def apply_transform(coords, R, mc, fc):
    return (np.asarray(coords) - mc) @ R + fc


def min_dist(a, b):
    if len(a) == 0 or len(b) == 0:
        return float("inf")
    return float(cKDTree(np.asarray(b)).query(np.asarray(a))[0].min())


def contact_count(binder_xyz, target_atoms):
    if len(target_atoms) == 0:
        return 0
    d = cKDTree(target_atoms.coord).query(binder_xyz)[0]
    return int((d < CONTACT_A).sum())


def ring_frame(target):
    ca = target[target.atom_name == "CA"]
    centre = ca.coord.mean(axis=0)
    xyz = ca.coord - centre

    u, s, vt = np.linalg.svd(
        xyz - xyz.mean(axis=0),
        full_matrices=False
    )

    axis = vt[2]
    axis = axis / np.linalg.norm(axis)

    return centre, axis


def rotation(axis, theta):
    a = axis / np.linalg.norm(axis)

    K = np.array([
        [0, -a[2], a[1]],
        [a[2], 0, -a[0]],
        [-a[1], a[0], 0],
    ])

    return (
        np.eye(3)
        + np.sin(theta) * K
        + (1 - np.cos(theta)) * (K @ K)
    )


def rotated_copy(xyz, centre, axis, theta):
    R = rotation(axis, theta)
    return (xyz - centre) @ R.T + centre


def build_alignment(st, meta, template):
    imap = meta.get("diffused_index_map", {})

    if not imap:
        raise ValueError("JSON has no diffused_index_map")

    dst_counts = Counter()

    for dst in imap.values():
        chain, resid = parse_label(dst)
        dst_counts[chain] += 1

    target_chain = dst_counts.most_common(1)[0][0]

    out_target = st[st.chain_id == target_chain]

    if len(out_target) == 0:
        raise ValueError(
            f"mapped target chain {target_chain} absent from CIF"
        )

    out_ca = {
        int(resid): xyz
        for resid, xyz in zip(
            out_target[out_target.atom_name == "CA"].res_id,
            out_target[out_target.atom_name == "CA"].coord,
        )
    }

    template_ca = {
        (str(chain), int(resid)): xyz
        for chain, resid, xyz in zip(
            template[template.atom_name == "CA"].chain_id,
            template[template.atom_name == "CA"].res_id,
            template[template.atom_name == "CA"].coord,
        )
    }

    moving = []
    fixed = []

    for src, dst in imap.items():
        src_chain, src_res = parse_label(src)
        dst_chain, dst_res = parse_label(dst)

        if dst_chain != target_chain:
            continue

        key = (src_chain, src_res)

        if key in template_ca and dst_res in out_ca:
            moving.append(out_ca[dst_res])
            fixed.append(template_ca[key])

    if len(moving) < 3:
        raise ValueError(
            f"only {len(moving)} CA matches available for alignment"
        )

    R, mc, fc, rmsd = kabsch(
        np.asarray(moving),
        np.asarray(fixed),
    )

    return target_chain, R, mc, fc, rmsd, len(moving)


def metric(meta, key):
    return meta.get("metrics", {}).get(key)


def duo_result(
    cif,
    meta,
    st,
    template,
    ring,
    centre,
    axis,
    target_chain,
    binder_chains,
    R,
    mc,
    fc,
    align_rmsd,
    n_match,
):
    if len(binder_chains) != 1:
        raise ValueError(
            f"duo expected 1 binder chain; found {binder_chains}"
        )

    binder_chain = binder_chains[0]
    binder = st[st.chain_id == binder_chain]
    xyz = apply_transform(binder.coord, R, mc, fc)

    nres = len(set(int(x) for x in binder.res_id))

    target_A = template[template.chain_id == "A"]
    target_B = template[template.chain_id == "B"]

    if len(target_A) == 0 or len(target_B) == 0:
        raise ValueError("duo target does not contain chains A and B")

    contacts_A = contact_count(xyz, target_A)
    contacts_B = contact_count(xyz, target_B)

    bridges = (
        contacts_A > MIN_CONTACTS
        and contacts_B > MIN_CONTACTS
    )

    local_target_min = min_dist(xyz, template.coord)
    ring_min = min_dist(xyz, ring.coord)

    nxt = rotated_copy(
        xyz, centre, axis, 2 * np.pi / 11
    )
    prv = rotated_copy(
        xyz, centre, axis, -2 * np.pi / 11
    )

    d_next = min_dist(xyz, nxt)
    d_prev = min_dist(xyz, prv)
    nbr_min = min(d_next, d_prev)

    target_clash = ring_min < CLASH_A
    symmetry_clash = nbr_min < MIN_NBR_A
    alignment_ok = align_rmsd < MAX_TARGET_CA_RMSD_A

    tight = (
        ring_min < COMFORT_A
        or nbr_min < COMFORT_A
    )

    passed = (
        alignment_ok
        and not target_clash
        and not symmetry_clash
        and bridges
    )

    return {
        "mode": "duo",
        "binder_chain": binder_chain,
        "binder_length": nres,
        "alignment_CA_count": n_match,
        "target_alignment_CA_RMSD_A": round(align_rmsd, 4),
        "alignment_ok": alignment_ok,
        "contacts_A": contacts_A,
        "contacts_B": contacts_B,
        "bridges_both_subunits": bridges,
        "local_target_min_dist_A": round(local_target_min, 3),
        "binder_ring_min_dist_A": round(ring_min, 3),
        "next_binder_min_dist_A": round(d_next, 3),
        "previous_binder_min_dist_A": round(d_prev, 3),
        "nbr_min_dist_A": round(nbr_min, 3),
        "target_clash": target_clash,
        "symmetry_clash": symmetry_clash,
        "tight_but_ok": tight,
        "rfd3_sidechain_clashes": metric(
            meta,
            "n_clashing.interresidue_clashes_w_sidechain",
        ),
        "rfd3_backbone_clashes": metric(
            meta,
            "n_clashing.interresidue_clashes_w_backbone",
        ),
        "PASS": passed,
    }


def trio_result(
    cif,
    meta,
    st,
    template,
    ring,
    centre,
    axis,
    target_chain,
    binder_chains,
    R,
    mc,
    fc,
    align_rmsd,
    n_match,
):
    if len(binder_chains) != 2:
        raise ValueError(
            f"trio expected 2 binder chains; found {binder_chains}"
        )

    for c in ("K", "A", "B"):
        if len(template[template.chain_id == c]) == 0:
            raise ValueError(
                f"trio target missing chain {c}"
            )

    xyz = {}
    lengths = {}

    for b in binder_chains:
        arr = st[st.chain_id == b]
        xyz[b] = apply_transform(arr.coord, R, mc, fc)
        lengths[b] = len(set(int(x) for x in arr.res_id))

    b0, b1 = binder_chains

    pair_min = min_dist(xyz[b0], xyz[b1])

    contacts = {}

    for b in binder_chains:
        for c in ("K", "A", "B"):
            sub = template[template.chain_id == c]
            contacts[(b, c)] = contact_count(xyz[b], sub)

    # Two possible assignments:
    # binder 1 = K/A and binder 2 = A/B, or vice versa
    assignments = [
        {
            b0: ("K", "A"),
            b1: ("A", "B"),
        },
        {
            b0: ("A", "B"),
            b1: ("K", "A"),
        },
    ]

    def assignment_score(assign):
        total = 0
        for b, pair in assign.items():
            total += contacts[(b, pair[0])]
            total += contacts[(b, pair[1])]
        return total

    assignment = max(
        assignments,
        key=assignment_score,
    )

    bridge_each = {}

    for b, pair in assignment.items():
        bridge_each[b] = (
            contacts[(b, pair[0])] > MIN_CONTACTS
            and contacts[(b, pair[1])] > MIN_CONTACTS
        )

    both_bridge = all(bridge_each.values())

    ring_min = {
        b: min_dist(xyz[b], ring.coord)
        for b in binder_chains
    }

    target_clash = any(
        ring_min[b] < CLASH_A
        for b in binder_chains
    )

    pair_clash = pair_min < MIN_NBR_A

    # Extra diagnostic:
    # would each backbone clear a +/-32.7 degree copy of itself?
    self_sym = {}

    for b in binder_chains:
        nxt = rotated_copy(
            xyz[b], centre, axis, 2 * np.pi / 11
        )
        prv = rotated_copy(
            xyz[b], centre, axis, -2 * np.pi / 11
        )

        self_sym[b] = min(
            min_dist(xyz[b], nxt),
            min_dist(xyz[b], prv),
        )

    alignment_ok = align_rmsd < MAX_TARGET_CA_RMSD_A

    tight = (
        pair_min < COMFORT_A
        or any(
            ring_min[b] < COMFORT_A
            for b in binder_chains
        )
    )

    passed = (
        alignment_ok
        and not target_clash
        and not pair_clash
        and both_bridge
    )

    out = {
        "mode": "trio_pair",
        "binder_chains": ",".join(binder_chains),
        "binder_length_1": lengths[b0],
        "binder_length_2": lengths[b1],
        "alignment_CA_count": n_match,
        "target_alignment_CA_RMSD_A": round(align_rmsd, 4),
        "alignment_ok": alignment_ok,
        "binder_binder_min_dist_A": round(pair_min, 3),
        "binder_binder_clash": pair_clash,
        "target_clash": target_clash,
        "tight_but_ok": tight,
        "assignment": (
            f"{b0}->{assignment[b0][0]}/{assignment[b0][1]};"
            f"{b1}->{assignment[b1][0]}/{assignment[b1][1]}"
        ),
        "both_bridge": both_bridge,
        "PASS": passed,
    }

    for b in binder_chains:
        out[f"{b}_ring_min_dist_A"] = round(ring_min[b], 3)
        out[f"{b}_self_symmetry_min_dist_A"] = round(
            self_sym[b], 3
        )
        out[f"{b}_self_symmetry_clash"] = (
            self_sym[b] < MIN_NBR_A
        )
        out[f"{b}_bridges_assigned_pair"] = bridge_each[b]

        for c in ("K", "A", "B"):
            out[f"contacts_{b}_{c}"] = contacts[(b, c)]

    return out


def check_one(cif, ring, centre, axis):
    cif = Path(cif)
    js = sidecar_for(cif)

    base = {
        "arm": cif.parent.name,
        "design": str(cif),
        "json": str(js),
    }

    try:
        if not js.exists():
            raise FileNotFoundError(
                f"missing sidecar JSON: {js}"
            )

        with open(js) as fh:
            meta = json.load(fh)

        spec = meta.get("specification", {})
        target_path = spec.get("input")

        if not target_path:
            raise ValueError(
                "JSON specification has no input target path"
            )

        target_path = Path(target_path)

        if not target_path.exists():
            raise FileNotFoundError(
                f"target input not found: {target_path}"
            )

        st = load_structure(cif)
        template = load_structure(target_path)

        (
            target_chain,
            R,
            mc,
            fc,
            align_rmsd,
            n_match,
        ) = build_alignment(
            st,
            meta,
            template,
        )

        chains = sorted(set(str(x) for x in st.chain_id))

        binder_chains = [
            c
            for c in chains
            if c != target_chain
        ]

        base["rfd3_target_chain"] = target_chain

        if len(binder_chains) == 1:
            result = duo_result(
                cif,
                meta,
                st,
                template,
                ring,
                centre,
                axis,
                target_chain,
                binder_chains,
                R,
                mc,
                fc,
                align_rmsd,
                n_match,
            )

        elif len(binder_chains) == 2:
            result = trio_result(
                cif,
                meta,
                st,
                template,
                ring,
                centre,
                axis,
                target_chain,
                binder_chains,
                R,
                mc,
                fc,
                align_rmsd,
                n_match,
            )

        else:
            raise ValueError(
                f"unsupported binder chains: {binder_chains}"
            )

        base.update(result)
        return base

    except Exception as exc:
        base["PASS"] = False
        base["error"] = f"{type(exc).__name__}: {exc}"
        return base


def write_csv(rows, path):
    keys = sorted(
        set().union(
            *(row.keys() for row in rows)
        )
    )

    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=keys,
        )
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--root",
        required=True,
    )
    ap.add_argument(
        "--ring",
        required=True,
    )
    ap.add_argument(
        "--json-out",
        default="out/rfd3_T20_screen/filter_results.json",
    )
    ap.add_argument(
        "--csv-out",
        default="out/rfd3_T20_screen/filter_results.csv",
    )

    args = ap.parse_args()

    root = Path(args.root)
    ring = load_structure(args.ring)
    centre, axis = ring_frame(ring)

    files = sorted(root.rglob("*.cif.gz"))

    rows = [
        check_one(
            cif,
            ring,
            centre,
            axis,
        )
        for cif in files
    ]

    with open(args.json_out, "w") as fh:
        json.dump(
            rows,
            fh,
            indent=2,
        )

    write_csv(
        rows,
        args.csv_out,
    )

    print("=== RFD3 C11 FILTER ===")
    print("designs checked:", len(rows))
    print()

    arms = sorted(set(r["arm"] for r in rows))

    for arm in arms:
        rr = [
            r
            for r in rows
            if r["arm"] == arm
        ]

        passed = sum(
            r.get("PASS") is True
            for r in rr
        )

        errors = sum(
            "error" in r
            for r in rr
        )

        print(
            f"{arm:22s} "
            f"{passed}/{len(rr)} PASS "
            f"errors={errors}"
        )

    print()
    print(
        "TOTAL PASS:",
        sum(r.get("PASS") is True for r in rows),
        "/",
        len(rows),
    )
    print("JSON:", args.json_out)
    print("CSV: ", args.csv_out)


if __name__ == "__main__":
    main()
