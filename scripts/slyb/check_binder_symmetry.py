#!/usr/bin/env python
"""Symmetry-expansion clash check for staple binders on the SlyB C11 ring.
CPU only. Usage: python check_binder_symmetry.py <design.pdb> [--order 11]

Takes a 3-chain complex (target A, target B, binder), places 11 symmetry
copies of the binder around the ring, and reports whether adjacent copies
clash. A design that scores well as a 1:1 binder but clashes with its own
neighbour is unusable -- this is the check that catches it.
"""
import sys, json, numpy as np
import biotite.structure as struc, biotite.structure.io.pdb as pdbio
from scipy.spatial import cKDTree

# Thresholds CALIBRATED against the native 7OJG interface, not assumed.
# The native A/B soluble-domain interface closes to 2.26 A with 13 heavy-atom
# pairs below 3.4 A -- a 3.4 A floor would reject the real protein. 2.2 A is
# set just below the native minimum: it still catches genuine interpenetration
# while admitting the close packing that real interfaces exhibit.
CLASH_A = 2.2          # binder-into-target hard clash floor
MIN_NBR_A = 2.2        # binder copy i vs copy i+1 (same rationale)
COMFORT_A = 3.4        # advisory: below this, flag as tight but not rejected

def load(p):
    """Read PDB, mmCIF or gzipped mmCIF.

    RFdiffusion3 writes <name>.cif.gz, so PDB-only parsing silently fails on
    every production design. Gzipped CIF is decompressed to a text handle
    rather than a temp file.
    """
    low = p.lower()
    if low.endswith((".cif", ".cif.gz", ".mmcif", ".mmcif.gz")):
        import biotite.structure.io.pdbx as pdbx
        if low.endswith(".gz"):
            import gzip, io
            with gzip.open(p, "rt") as fh:
                f = pdbx.CIFFile.read(io.StringIO(fh.read()))
        else:
            f = pdbx.CIFFile.read(p)
        a = pdbx.get_structure(f, model=1)
    else:
        a = pdbio.PDBFile.read(p).get_structure(model=1)
    return a[(a.element != "H") & struc.filter_amino_acids(a)]

def ring_frame(target):
    """Recover axis/centre from the 11-fold target itself."""
    ca = target[target.atom_name == "CA"]
    centre = ca.coord.mean(0)
    xyz = ca.coord - centre
    # pore axis = smallest-inertia direction of a flat ring is the normal
    u, s, vt = np.linalg.svd(xyz - xyz.mean(0), full_matrices=False)
    axis = vt[2] / np.linalg.norm(vt[2])
    return centre, axis

def rot(axis, th):
    a = axis / np.linalg.norm(axis)
    K = np.array([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])
    return np.eye(3) + np.sin(th)*K + (1-np.cos(th))*(K@K)

def check_pair(path, target_chains=("K","A","B")):
    """Two-binder trio test: both binders explicitly present in one file.
    Reports binder-binder distance and per-target-chain contacts, so an
    asymmetric or clashing pair is caught without symmetry expansion."""
    st = load(path)
    bch = [c for c in sorted(set(st.chain_id)) if c not in target_chains]
    if len(bch) != 2:
        return dict(design=path, error="expected exactly 2 binder chains, found %s" % bch)
    b1 = st[st.chain_id == bch[0]]; b2 = st[st.chain_id == bch[1]]
    tgt = st[np.isin(st.chain_id, list(target_chains))]
    d_bb = float(cKDTree(b2.coord).query(b1.coord)[0].min())
    out = dict(design=path, binder_chains=bch, binder_binder_min_A=round(d_bb, 2),
               binder_binder_clash=bool(d_bb < MIN_NBR_A))
    for lab, b in ((bch[0], b1), (bch[1], b2)):
        d,_ = cKDTree(tgt.coord).query(b.coord)
        out["min_to_target_%s" % lab] = round(float(d.min()), 2)
        for c in target_chains:
            sub = st[st.chain_id == c]
            if len(sub) == 0: continue
            dd,_ = cKDTree(sub.coord).query(b.coord)
            out["contacts_%s_%s" % (lab, c)] = int((dd < 5.0).sum())
    bridging = all(
        sum(1 for c in target_chains if out.get("contacts_%s_%s" % (lab, c), 0) > 20) >= 2
        for lab in bch)
    clash = any(out["min_to_target_%s" % l] < CLASH_A for l in bch)
    out["both_bridge"] = bool(bridging)
    out["tight_but_ok"] = bool(d_bb < COMFORT_A or
                               any(out["min_to_target_%s" % l] < COMFORT_A for l in bch))
    out["PASS"] = bool(not out["binder_binder_clash"] and not clash and bridging)
    return out


def check(path, order=11, binder_chain=None, ref_ring=None):
    st = load(path)
    chains = sorted(set(st.chain_id))
    if binder_chain is None:
        # binder = the chain that is not A/B
        cand = [c for c in chains if c not in ("A","B")]
        if not cand:
            return dict(design=path, error="no binder chain found (expected a chain besides A/B)")
        binder_chain = cand[0]
    binder = st[st.chain_id == binder_chain]
    target = st[st.chain_id != binder_chain]
    if len(binder) == 0 or len(target) == 0:
        return dict(design=path, error="empty binder or target selection")

    if ref_ring is not None:
        centre, axis = ring_frame(load(ref_ring))
    else:
        centre, axis = ring_frame(target)

    R = rot(axis, 2*np.pi/order)
    nxt = (binder.coord - centre) @ R.T + centre
    prv = (binder.coord - centre) @ rot(axis, -2*np.pi/order).T + centre
    d_next = float(cKDTree(nxt).query(binder.coord)[0].min())
    d_prev = float(cKDTree(prv).query(binder.coord)[0].min())

    # binder vs target clash
    d_tgt = float(cKDTree(target.coord).query(binder.coord)[0].min())
    # interface size: binder heavy atoms within 5 A of each target chain
    res = {}
    for c in sorted(set(target.chain_id)):
        sub = target[target.chain_id == c]
        d,_ = cKDTree(sub.coord).query(binder.coord)
        res["contacts_%s" % c] = int((d < 5.0).sum())

    bridges_both = res.get("contacts_A",0) > 20 and res.get("contacts_B",0) > 20
    out = dict(design=path, binder_chain=binder_chain, n_binder_atoms=len(binder),
               nbr_min_dist_A=round(min(d_next,d_prev),2),
               binder_target_min_dist_A=round(d_tgt,2),
               symmetry_clash=bool(min(d_next,d_prev) < MIN_NBR_A),
               target_clash=bool(d_tgt < CLASH_A),
               bridges_both_subunits=bool(bridges_both), **res)
    out["tight_but_ok"] = bool(min(d_next, d_prev) < COMFORT_A or d_tgt < COMFORT_A)
    out["PASS"] = bool(not out["symmetry_clash"] and not out["target_clash"]
                       and out["bridges_both_subunits"])
    return out

if __name__ == "__main__":
    argv = sys.argv[1:]
    order, ref, args, i = 11, None, [], 0
    while i < len(argv):
        if argv[i] == "--order":
            order = int(argv[i+1]); i += 2
        elif argv[i] == "--ref":
            ref = argv[i+1]; i += 2
        elif argv[i].startswith("--"):
            i += 1
        else:
            args.append(argv[i]); i += 1
    if not args:
        sys.exit("usage: check_binder_symmetry.py <design.pdb> [...] "
                 "[--order N] [--ref ring.pdb] [--pair]")
    if "--pair" in argv:
        rows = [check_pair(p) for p in args]
    else:
        rows = [check(p, order=order, ref_ring=ref) for p in args]
    print(json.dumps(rows, indent=2))
