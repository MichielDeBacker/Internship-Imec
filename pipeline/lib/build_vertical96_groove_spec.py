#!/usr/bin/env python3

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import biotite.structure.io.pdb as pdbio
import biotite.structure.io.pdbx as pdbx


ALLOWED_96 = set(range(18, 61)) | set(range(103, 156))
BACKBONE = {"N", "CA", "C", "O", "OXT"}


def load(path):
    path = Path(path)

    if path.suffix.lower() == ".pdb":
        f = pdbio.PDBFile.read(path)
        return f.get_structure(model=1)

    f = pdbx.CIFFile.read(path)
    return pdbx.get_structure(f, model=1)


def write_pdb(a, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    f = pdbio.PDBFile()
    f.set_structure(a)
    f.write(path)


def chain_ids(a):
    result = []

    for cid in a.chain_id:
        cid = str(cid)

        if cid not in result:
            result.append(cid)

    return result


def get_chain(a, cid):
    mask = np.asarray(
        [str(x) == str(cid) for x in a.chain_id],
        dtype=bool,
    )
    return a[mask]


def crop96(a):
    mask = np.asarray(
        [int(x) in ALLOWED_96 for x in a.res_id],
        dtype=bool,
    )
    return a[mask]


def ca(a):
    mask = np.asarray(
        [str(x) == "CA" for x in a.atom_name],
        dtype=bool,
    )
    return np.asarray(a.coord[mask], dtype=float)


def kabsch(mobile, fixed):
    mobile = np.asarray(mobile, dtype=float)
    fixed = np.asarray(fixed, dtype=float)

    if mobile.shape != fixed.shape:
        raise ValueError(
            f"Kabsch mismatch {mobile.shape} vs {fixed.shape}"
        )

    if len(mobile) < 3:
        raise ValueError("Need >=3 coordinates")

    cm = mobile.mean(axis=0)
    cf = fixed.mean(axis=0)

    x = mobile - cm
    y = fixed - cf

    u, _, vt = np.linalg.svd(x.T @ y)

    r = vt.T @ u.T

    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T

    t = cf - cm @ r.T

    fit = mobile @ r.T + t

    rmsd = float(
        np.sqrt(
            np.mean(
                np.sum(
                    (fit - fixed) ** 2,
                    axis=1,
                )
            )
        )
    )

    return r, t, rmsd


def transform(a, r, t):
    b = a.copy()

    b.coord = (
        np.asarray(a.coord, dtype=float) @ r.T
        + t
    )

    return b


def ring_geometry(ring):
    ids = chain_ids(ring)

    coms = np.asarray(
        [
            ca(get_chain(ring, cid)).mean(axis=0)
            for cid in ids
        ]
    )

    centre = coms.mean(axis=0)
    centered = coms - centre

    vals, vecs = np.linalg.eigh(
        centered.T @ centered
    )

    axis_i = int(np.argmin(vals))

    axis = vecs[:, axis_i]
    axis /= np.linalg.norm(axis)

    plane = [
        i
        for i in np.argsort(vals)[::-1]
        if i != axis_i
    ]

    e1 = vecs[:, plane[0]]
    e2 = vecs[:, plane[1]]

    rows = []

    for cid, xyz in zip(ids, centered):
        angle = math.atan2(
            float(xyz @ e2),
            float(xyz @ e1),
        )
        rows.append((angle, cid))

    rows.sort()

    return (
        [x[1] for x in rows],
        centre,
        axis,
    )


def match_penta_to_ring(penta, ring):
    mobile = []
    fixed = []

    ring_ids = set(chain_ids(ring))

    for cid in chain_ids(penta):
        if cid not in ring_ids:
            continue

        pc = get_chain(penta, cid)
        rc = get_chain(ring, cid)

        ref = {
            int(rc.res_id[i]):
            np.asarray(rc.coord[i], dtype=float)
            for i in range(len(rc))
            if str(rc.atom_name[i]) == "CA"
        }

        for i in range(len(pc)):
            if str(pc.atom_name[i]) != "CA":
                continue

            rid = int(pc.res_id[i])

            if rid in ref:
                mobile.append(pc.coord[i])
                fixed.append(ref[rid])

    if len(mobile) < 100:
        raise ValueError(
            f"Only {len(mobile)} penta/native CA matches"
        )

    return kabsch(
        np.asarray(mobile),
        np.asarray(fixed),
    )


def binder_contact_count(
    binder,
    target,
    cutoff=5.0,
):
    tree = cKDTree(
        np.asarray(target.coord, dtype=float)
    )

    distances, _ = tree.query(
        np.asarray(binder.coord, dtype=float),
        k=1,
    )

    binder_residues = set()

    for i, d in enumerate(distances):
        if d <= cutoff:
            binder_residues.add(
                int(binder.res_id[i])
            )

    return len(binder_residues)


def hotspot_rows(
    target,
    binder,
    axis,
    centre,
    max_distance=6.0,
):
    """
    Candidate hotspots:
      - only historical96 residues
      - side-chain heavy atoms
      - physically close to mapped historical binder
    """

    binder_tree = cKDTree(
        np.asarray(binder.coord, dtype=float)
    )

    by_residue = {}

    for i in range(len(target)):
        rid = int(target.res_id[i])

        if rid not in ALLOWED_96:
            continue

        by_residue.setdefault(rid, []).append(i)

    rows = []

    for rid, idxs in by_residue.items():

        ca_idxs = [
            i
            for i in idxs
            if str(target.atom_name[i]) == "CA"
        ]

        if not ca_idxs:
            continue

        sidechain = [
            i
            for i in idxs
            if str(target.atom_name[i]) not in BACKBONE
            and not str(target.atom_name[i]).upper().startswith("H")
        ]

        # Do not use backbone-only/Gly positions as hotspots.
        if not sidechain:
            continue

        xyz = np.asarray(
            target.coord[sidechain],
            dtype=float,
        )

        distances, _ = binder_tree.query(
            xyz,
            k=1,
        )

        j = int(np.argmin(distances))
        atom_i = sidechain[j]
        distance = float(distances[j])

        if distance > max_distance:
            continue

        projection = float(
            (
                np.asarray(
                    target.coord[ca_idxs[0]],
                    dtype=float,
                )
                - centre
            )
            @ axis
        )

        rows.append(
            {
                "resid": rid,
                "resname":
                    str(target.res_name[atom_i]),
                "atom":
                    str(target.atom_name[atom_i]),
                "distance_A": distance,
                "axis_projection_A":
                    projection,
            }
        )

    rows.sort(
        key=lambda x: x["distance_A"]
    )

    return rows


def choose_two_hotspots(rows):
    """
    Pick two true interface residues with useful
    separation along the C11 symmetry axis.
    """

    if len(rows) < 2:
        raise ValueError(
            "Fewer than two groove-facing "
            "side-chain hotspot candidates"
        )

    # Restrict to nearest interface residues first.
    near = rows[:12]

    best = None

    for a, b in itertools.combinations(near, 2):

        axial = abs(
            a["axis_projection_A"]
            - b["axis_projection_A"]
        )

        # Prefer vertical span, but penalize
        # unnecessarily distant interface atoms.
        score = (
            axial
            - 0.5 * (
                a["distance_A"]
                + b["distance_A"]
            )
        )

        if best is None or score > best[0]:
            best = (
                score,
                axial,
                a,
                b,
            )

    _, axial, a, b = best

    # Need meaningful axial separation.
    if axial < 5.0:
        # Fall back to two closest true interface
        # residues, but retain explicit warning.
        a, b = near[0], near[1]
        axial = abs(
            a["axis_projection_A"]
            - b["axis_projection_A"]
        )

    return [a, b], axial


def principal_axis_angle(binder, ring_axis):
    xyz = ca(binder)

    xyz = xyz - xyz.mean(axis=0)

    _, _, vt = np.linalg.svd(
        xyz,
        full_matrices=False,
    )

    binder_axis = vt[0]
    binder_axis /= np.linalg.norm(
        binder_axis
    )

    dot = abs(
        float(
            binder_axis
            @ ring_axis
        )
    )

    dot = max(-1.0, min(1.0, dot))

    return float(
        np.degrees(
            np.arccos(dot)
        )
    )


def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--ring",
        required=True,
    )

    ap.add_argument(
        "--penta",
        required=True,
    )

    ap.add_argument(
        "--parent",
        required=True,
    )

    ap.add_argument(
        "--out",
        required=True,
    )

    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    ring = load(args.ring)
    penta = load(args.penta)
    parent = load(args.parent)

    order, centre, axis = (
        ring_geometry(ring)
    )

    if len(order) != 11:
        raise SystemExit(
            f"Expected C11 ring, found "
            f"{len(order)} chains: {order}"
        )

    if "A" not in order or "B" not in order:
        raise SystemExit(
            f"Native ring lacks A/B: {order}"
        )

    ia = order.index("A")
    ib = order.index("B")

    if order[
        (ia + 1) % len(order)
    ] == "B":
        ref_first = "A"
        ref_second = "B"

    elif order[
        (ib + 1) % len(order)
    ] == "A":
        ref_first = "B"
        ref_second = "A"

    else:
        raise SystemExit(
            "A/B are not adjacent in "
            f"native ring order: {order}"
        )

    # ------------------------------------------
    # Recover design2 historical binder pose
    # ------------------------------------------

    pids = chain_ids(parent)

    if len(pids) != 2:
        raise SystemExit(
            f"design2 should have 2 chains, "
            f"found {pids}"
        )

    lengths = {
        cid: len(
            ca(
                get_chain(
                    parent,
                    cid,
                )
            )
        )
        for cid in pids
    }

    binder_id = min(
        pids,
        key=lambda x: lengths[x],
    )

    target_id = max(
        pids,
        key=lambda x: lengths[x],
    )

    binder = get_chain(
        parent,
        binder_id,
    )

    historical_target = get_chain(
        parent,
        target_id,
    )

    if len(ca(historical_target)) != 480:
        raise SystemExit(
            "Expected historical merged target "
            f"to contain 480 CA atoms, found "
            f"{len(ca(historical_target))}"
        )

    # Exact reconstruction of penta concatenation
    best = None

    for perm in itertools.permutations(
        chain_ids(penta)
    ):
        reference = np.concatenate(
            [
                ca(
                    get_chain(
                        penta,
                        cid,
                    )
                )
                for cid in perm
            ],
            axis=0,
        )

        if reference.shape != ca(
            historical_target
        ).shape:
            continue

        r, t, rmsd = kabsch(
            ca(historical_target),
            reference,
        )

        if best is None or rmsd < best[0]:
            best = (
                rmsd,
                perm,
                r,
                t,
            )

    if best is None:
        raise SystemExit(
            "Could not reconstruct "
            "design2 -> penta mapping"
        )

    parent_rmsd, penta_order, r1, t1 = best

    if parent_rmsd > 2.0:
        raise SystemExit(
            f"design2 -> penta RMSD "
            f"{parent_rmsd:.3f} A > 2 A"
        )

    binder_penta = transform(
        binder,
        r1,
        t1,
    )

    r2, t2, penta_rmsd = (
        match_penta_to_ring(
            penta,
            ring,
        )
    )

    if penta_rmsd > 2.0:
        raise SystemExit(
            f"penta -> native C11 RMSD "
            f"{penta_rmsd:.3f} A > 2 A"
        )

    binder_ring = transform(
        binder_penta,
        r2,
        t2,
    )

    # ------------------------------------------
    # Determine historical adjacent groove
    # ------------------------------------------

    candidates = []

    for i, c1 in enumerate(order):

        c2 = order[
            (i + 1) % len(order)
        ]

        n1 = binder_contact_count(
            binder_ring,
            crop96(
                get_chain(
                    ring,
                    c1,
                )
            ),
        )

        n2 = binder_contact_count(
            binder_ring,
            crop96(
                get_chain(
                    ring,
                    c2,
                )
            ),
        )

        candidates.append(
            (
                min(n1, n2),
                n1 + n2,
                c1,
                c2,
                n1,
                n2,
            )
        )

    candidates.sort(
        reverse=True
    )

    (
        _,
        _,
        hist1,
        hist2,
        hist_contacts1,
        hist_contacts2,
    ) = candidates[0]

    if min(
        hist_contacts1,
        hist_contacts2,
    ) < 3:
        raise SystemExit(
            "Historical binder does not make "
            ">=3 binder-residue contacts to "
            "both cropped96 partners"
        )

    # ------------------------------------------
    # Map historical groove -> canonical A/B
    # ------------------------------------------

    src_pair = np.concatenate(
        [
            ca(
                get_chain(
                    ring,
                    hist1,
                )
            ),
            ca(
                get_chain(
                    ring,
                    hist2,
                )
            ),
        ],
        axis=0,
    )

    dst_pair = np.concatenate(
        [
            ca(
                get_chain(
                    ring,
                    ref_first,
                )
            ),
            ca(
                get_chain(
                    ring,
                    ref_second,
                )
            ),
        ],
        axis=0,
    )

    rg, tg, groove_rmsd = kabsch(
        src_pair,
        dst_pair,
    )

    if groove_rmsd > 1.0:
        raise SystemExit(
            f"historical -> canonical groove "
            f"RMSD {groove_rmsd:.3f} A > 1 A"
        )

    binder_canonical = transform(
        binder_ring,
        rg,
        tg,
    )

    target1 = crop96(
        get_chain(
            ring,
            ref_first,
        )
    )

    target2 = crop96(
        get_chain(
            ring,
            ref_second,
        )
    )

    if len(ca(target1)) != 96:
        raise SystemExit(
            f"{ref_first}: expected 96 CA, "
            f"found {len(ca(target1))}"
        )

    if len(ca(target2)) != 96:
        raise SystemExit(
            f"{ref_second}: expected 96 CA, "
            f"found {len(ca(target2))}"
        )

    # ------------------------------------------
    # Select genuine groove-facing hotspots
    # ------------------------------------------

    rows1 = hotspot_rows(
        target1,
        binder_canonical,
        axis,
        centre,
    )

    rows2 = hotspot_rows(
        target2,
        binder_canonical,
        axis,
        centre,
    )

    h1, span1 = choose_two_hotspots(
        rows1
    )

    h2, span2 = choose_two_hotspots(
        rows2
    )

    hotspots = {}

    for row in h1:
        hotspots[
            f"{ref_first}{row['resid']}"
        ] = row["atom"]

    for row in h2:
        hotspots[
            f"{ref_second}{row['resid']}"
        ] = row["atom"]

    if len(hotspots) != 4:
        raise SystemExit(
            f"Expected 4 distinct hotspots, "
            f"found {hotspots}"
        )

    # All hotspots must belong to historical96.
    for key in hotspots:
        rid = int(key[1:])

        if rid not in ALLOWED_96:
            raise SystemExit(
                f"Illegal hotspot outside "
                f"historical96: {key}"
            )

    canonical_contacts1 = (
        binder_contact_count(
            binder_canonical,
            target1,
        )
    )

    canonical_contacts2 = (
        binder_contact_count(
            binder_canonical,
            target2,
        )
    )

    if min(
        canonical_contacts1,
        canonical_contacts2,
    ) < 3:
        raise SystemExit(
            "Mapped canonical binder does not "
            "contact both 96-aa partners"
        )

    vertical_angle = (
        principal_axis_angle(
            binder_canonical,
            axis,
        )
    )

    # ------------------------------------------
    # Write actual 96 + 96 target
    # ------------------------------------------

    target96 = target1 + target2

    target_path = (
        out / "target_AB_historical96.pdb"
    )

    write_pdb(
        target96,
        target_path,
    )

    # ------------------------------------------
    # Local visual preview
    # ------------------------------------------

    preview_binder = (
        binder_canonical.copy()
    )

    preview_binder.chain_id[:] = "X"

    write_pdb(
        target96 + preview_binder,
        out / "preview_AB96_plus_design2_binder.pdb",
    )

    # ------------------------------------------
    # Full C11 + 11 mapped binder copies
    # ------------------------------------------

    full_preview = ring.copy()

    copy_ids = list(
        "0123456789Z"
    )

    canonical_pair = np.concatenate(
        [
            ca(
                get_chain(
                    ring,
                    ref_first,
                )
            ),
            ca(
                get_chain(
                    ring,
                    ref_second,
                )
            ),
        ],
        axis=0,
    )

    for i, c1 in enumerate(order):

        c2 = order[
            (i + 1) % len(order)
        ]

        pair = np.concatenate(
            [
                ca(
                    get_chain(
                        ring,
                        c1,
                    )
                ),
                ca(
                    get_chain(
                        ring,
                        c2,
                    )
                ),
            ],
            axis=0,
        )

        rr, tt, _ = kabsch(
            canonical_pair,
            pair,
        )

        bcopy = transform(
            binder_canonical,
            rr,
            tt,
        )

        bcopy.chain_id[:] = copy_ids[i]

        full_preview = (
            full_preview + bcopy
        )

    write_pdb(
        full_preview,
        out / "preview_full_C11_plus_11_binders.pdb",
    )

    # ------------------------------------------
    # RFD3 specs
    # ------------------------------------------

    contig_target = (
        f"{ref_first}18-60,"
        f"{ref_first}103-155,"
        f"/0,"
        f"{ref_second}18-60,"
        f"{ref_second}103-155"
    )

    specs = {}

    for name, length in [
        (
            "vertical96_46_54",
            "46-54",
        ),
        (
            "vertical96_54_62",
            "54-62",
        ),
        (
            "vertical96_62_72",
            "62-72",
        ),
    ]:

        specs[name] = {
            "dialect": 2,

            "input":
                str(
                    target_path.resolve()
                ),

            "contig":
                f"{length},/0,"
                f"{contig_target}",

            "infer_ori_strategy":
                "hotspots",

            "select_hotspots":
                hotspots,

            "is_non_loopy":
                True,

            "plddt_enhanced":
                True,
        }

    spec_path = (
        out / "vertical96_spec.json"
    )

    spec_path.write_text(
        json.dumps(
            specs,
            indent=2,
        )
        + "\n"
    )

    geometry = {
        "target_representation":
            "historical96",

        "target_residues":
            [
                "18-60",
                "103-155",
            ],

        "target_residues_per_chain":
            96,

        "ring_order":
            order,

        "ring_axis":
            axis.tolist(),

        "ring_centre":
            centre.tolist(),

        "historical_parent":
            "design2",

        "historical_binder_chain":
            binder_id,

        "historical_target_chain":
            target_id,

        "historical_groove":
            [
                hist1,
                hist2,
            ],

        "generation_groove":
            [
                ref_first,
                ref_second,
            ],

        "historical_contacts96":
            {
                hist1:
                    hist_contacts1,

                hist2:
                    hist_contacts2,
            },

        "canonical_contacts96":
            {
                ref_first:
                    canonical_contacts1,

                ref_second:
                    canonical_contacts2,
            },

        "hotspots":
            hotspots,

        "hotspot_details":
            {
                ref_first: h1,
                ref_second: h2,
            },

        "hotspot_axis_span_A":
            {
                ref_first:
                    span1,

                ref_second:
                    span2,
            },

        "historical_binder_vertical_angle_deg":
            vertical_angle,

        "design2_to_penta_RMSD_A":
            parent_rmsd,

        "penta_to_native_C11_RMSD_A":
            penta_rmsd,

        "historical_to_generation_groove_RMSD_A":
            groove_rmsd,

        "penta_chain_order":
            list(penta_order),
    }

    geo_path = (
        out / "vertical96_geometry.json"
    )

    geo_path.write_text(
        json.dumps(
            geometry,
            indent=2,
        )
        + "\n"
    )

    print()
    print(
        "=========================================="
    )
    print(
        "HISTORICAL96 VERTICAL GROOVE READY"
    )
    print(
        "=========================================="
    )

    print(
        "ring order:",
        order,
    )

    print(
        "historical groove:",
        hist1,
        hist2,
    )

    print(
        "generation groove:",
        ref_first,
        ref_second,
    )

    print(
        "target CA:",
        len(ca(target1)),
        "+",
        len(ca(target2)),
    )

    print(
        "contacts:",
        canonical_contacts1,
        canonical_contacts2,
    )

    print()
    print(
        "SELECTED GROOVE HOTSPOTS:"
    )

    for chain, rows in [
        (ref_first, h1),
        (ref_second, h2),
    ]:
        for row in rows:
            print(
                f"  {chain}"
                f"{row['resid']} "
                f"{row['resname']} "
                f"atom={row['atom']} "
                f"historical_distance="
                f"{row['distance_A']:.2f} A "
                f"axis="
                f"{row['axis_projection_A']:.2f} A"
            )

    print()
    print(
        "RFD3 hotspot dictionary:",
        hotspots,
    )

    print(
        "design2->penta RMSD:",
        round(
            parent_rmsd,
            4,
        ),
    )

    print(
        "penta->C11 RMSD:",
        round(
            penta_rmsd,
            4,
        ),
    )

    print(
        "groove symmetry RMSD:",
        round(
            groove_rmsd,
            4,
        ),
    )

    print(
        "historical binder vertical angle:",
        round(
            vertical_angle,
            2,
        ),
        "deg",
    )

    print()
    print(
        "spec:",
        spec_path,
    )

    print(
        "AB preview:",
        out
        / "preview_AB96_plus_design2_binder.pdb",
    )

    print(
        "C11 preview:",
        out
        / "preview_full_C11_plus_11_binders.pdb",
    )


if __name__ == "__main__":
    main()
