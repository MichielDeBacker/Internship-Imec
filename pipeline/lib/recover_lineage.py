#!/usr/bin/env python3

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

import biotite.structure.io.pdb as pdbio
import biotite.structure.io.pdbx as pdbx


AA3 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
    "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
    "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}


def load(path):
    path = Path(path)
    low = path.name.lower()

    if low.endswith((".cif.gz", ".mmcif.gz")):
        with gzip.open(path, "rt") as fh:
            obj = pdbx.CIFFile.read(
                io.StringIO(fh.read())
            )
        return pdbx.get_structure(obj, model=1)

    if low.endswith((".cif", ".mmcif")):
        obj = pdbx.CIFFile.read(path)
        return pdbx.get_structure(obj, model=1)

    if low.endswith(".pdb"):
        obj = pdbio.PDBFile.read(path)
        return obj.get_structure(model=1)

    raise ValueError(path)


def chain_ids(a):
    out = []

    for c in a.chain_id:
        c = str(c)
        if c not in out:
            out.append(c)

    return out


def chain(a, cid):
    mask = np.asarray(
        [str(x) == str(cid) for x in a.chain_id],
        dtype=bool,
    )
    return a[mask]


def chain_data(a, cid):
    x = chain(a, cid)

    coords = []
    seq = []
    resids = []
    seen = set()

    for i in range(len(x)):
        if str(x.atom_name[i]) != "CA":
            continue

        rid = int(x.res_id[i])

        if rid in seen:
            continue

        seen.add(rid)

        coords.append(
            np.asarray(x.coord[i], dtype=float)
        )

        resids.append(rid)

        seq.append(
            AA3.get(
                str(x.res_name[i]).upper(),
                "X",
            )
        )

    return {
        "chain": str(cid),
        "n_res": len(coords),
        "coords": np.asarray(coords, dtype=float),
        "seq": "".join(seq),
        "resids": resids,
    }


def structure_data(path):
    a = load(path)

    return [
        chain_data(a, cid)
        for cid in chain_ids(a)
    ]


def seq_hash(seq):
    return hashlib.sha1(
        seq.encode()
    ).hexdigest()[:12]


def seq_identity(a, b):
    if len(a) != len(b) or not a:
        return 0.0

    return sum(
        x == y
        for x, y in zip(a, b)
    ) / len(a)


def kabsch_rmsd(a, b):
    if a.shape != b.shape:
        return None

    if len(a) < 3:
        return None

    ac = a.mean(axis=0)
    bc = b.mean(axis=0)

    x = a - ac
    y = b - bc

    u, s, vt = np.linalg.svd(
        x.T @ y
    )

    r = vt.T @ u.T

    if np.linalg.det(r) < 0:
        vt[-1] *= -1
        r = vt.T @ u.T

    fit = x @ r.T

    return float(
        np.sqrt(
            np.mean(
                np.sum(
                    (fit - y) ** 2,
                    axis=1,
                )
            )
        )
    )


def sidecar_for(path):
    s = str(path)

    if s.endswith(".cif.gz"):
        return Path(s[:-7] + ".json")

    if s.endswith(".cif"):
        return Path(s[:-4] + ".json")

    return None


def load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except Exception:
        return None


def nested_value(obj, key):
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]

        for value in obj.values():
            hit = nested_value(
                value,
                key,
            )

            if hit is not None:
                return hit

    elif isinstance(obj, list):
        for value in obj:
            hit = nested_value(
                value,
                key,
            )

            if hit is not None:
                return hit

    return None


def write_tsv(path, rows, fields):
    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with path.open(
        "w",
        newline="",
    ) as f:
        w = csv.DictWriter(
            f,
            delimiter="\t",
            lineterminator="\n",
            fieldnames=fields,
        )

        w.writeheader()

        for row in rows:
            w.writerow(
                {
                    field: row.get(
                        field,
                        "",
                    )
                    for field in fields
                }
            )


def candidate_name(path):
    m = re.search(
        r"(d\d+_s\d+)",
        str(path),
    )

    return m.group(1) if m else ""


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--slyb",
        required=True,
    )

    ap.add_argument(
        "--repo",
        required=True,
    )

    ap.add_argument(
        "--out",
        required=True,
    )

    args = ap.parse_args()

    slyb = Path(args.slyb)
    repo = Path(args.repo)
    out = Path(args.out)

    out.mkdir(
        parents=True,
        exist_ok=True,
    )

    parents = {
        "design0":
            slyb / "sequence_design/inputs/design0.cif",
        "design1":
            slyb / "sequence_design/inputs/design1.cif",
        "design2":
            slyb / "sequence_design/inputs/design2.cif",
        "design6":
            slyb / "sequence_design/inputs/design6.cif",
    }

    candidates = {
        "d2_s0":
            slyb / "sequence_design/out/rf3_primary_50056967/d2_s0/d2_s0_model.cif",
        "d2_s1":
            slyb / "sequence_design/out/rf3_primary_50056967/d2_s1/d2_s1_model.cif",
        "d1_s3":
            slyb / "sequence_design/out/rf3_primary_50056967/d1_s3/d1_s3_model.cif",
        "d6_s6":
            slyb / "sequence_design/out/rf3_primary_50056967/d6_s6/d6_s6_model.cif",
        "d6_s3":
            slyb / "sequence_design/out/rf3_primary_50056967/d6_s3/d6_s3_model.cif",
        "d6_s7":
            slyb / "sequence_design/out/rf3_primary_50056967/d6_s7/d6_s7_model.cif",
    }

    family = {
        "d2_s0": "design2",
        "d2_s1": "design2",
        "d1_s3": "design1",
        "d6_s6": "design6",
        "d6_s3": "design6",
        "d6_s7": "design6",
    }

    # --------------------------------------------------------
    # Native target
    # --------------------------------------------------------

    ring = (
        slyb
        / "inputs"
        / "target_ring_C11.pdb"
    )

    ring_data = structure_data(ring)

    native_a = next(
        x
        for x in ring_data
        if x["chain"] == "A"
    )

    print(
        "native SlyB chain A:",
        native_a["n_res"],
        "residues",
    )

    # --------------------------------------------------------
    # Parent structures
    # --------------------------------------------------------

    parent_data = {}

    structure_rows = []

    for name, path in parents.items():
        data = structure_data(path)
        parent_data[name] = data

        for c in data:
            structure_rows.append(
                {
                    "group": "parent",
                    "name": name,
                    "path": str(path),
                    "chain": c["chain"],
                    "n_res": c["n_res"],
                    "seq_sha1": seq_hash(
                        c["seq"]
                    ),
                    "native_seq_identity": (
                        seq_identity(
                            c["seq"],
                            native_a["seq"],
                        )
                    ),
                }
            )

    # --------------------------------------------------------
    # RF3 monomers
    # --------------------------------------------------------

    candidate_data = {}

    for name, path in candidates.items():
        data = structure_data(path)
        candidate_data[name] = data

        for c in data:
            structure_rows.append(
                {
                    "group": "rf3_monomer",
                    "name": name,
                    "path": str(path),
                    "chain": c["chain"],
                    "n_res": c["n_res"],
                    "seq_sha1": seq_hash(
                        c["seq"]
                    ),
                    "native_seq_identity": (
                        seq_identity(
                            c["seq"],
                            native_a["seq"],
                        )
                    ),
                }
            )

    # --------------------------------------------------------
    # Parent <-> monomer structural matching
    # --------------------------------------------------------

    monomer_matches = []

    for candidate, fam in family.items():
        mono = candidate_data[candidate]

        if len(mono) != 1:
            continue

        mono = mono[0]

        for pc in parent_data[fam]:
            rmsd = kabsch_rmsd(
                mono["coords"],
                pc["coords"],
            )

            monomer_matches.append(
                {
                    "candidate": candidate,
                    "family": fam,
                    "parent_chain": pc["chain"],
                    "parent_n_res": pc["n_res"],
                    "monomer_n_res": mono["n_res"],
                    "sequence_identity": round(
                        seq_identity(
                            mono["seq"],
                            pc["seq"],
                        ),
                        6,
                    ),
                    "CA_RMSD_A": (
                        ""
                        if rmsd is None
                        else round(
                            rmsd,
                            5,
                        )
                    ),
                }
            )

    # --------------------------------------------------------
    # Classify parent representation
    # --------------------------------------------------------

    family_monomer_lengths = (
        defaultdict(set)
    )

    for candidate, fam in family.items():
        if candidate_data[candidate]:
            family_monomer_lengths[
                fam
            ].add(
                candidate_data[
                    candidate
                ][0]["n_res"]
            )

    parent_class = {}

    for fam, chains in parent_data.items():
        lengths = [
            c["n_res"]
            for c in chains
        ]

        known = (
            family_monomer_lengths.get(
                fam,
                set(),
            )
        )

        if (
            len(chains) == 2
            and known
            and all(
                n in known
                for n in lengths
            )
        ):
            classification = (
                "two_binder_parent"
            )

        elif (
            len(chains) == 2
            and 138 in lengths
        ):
            classification = (
                "target_plus_binder"
            )

        else:
            classification = (
                "ambiguous_two_chain"
            )

        parent_class[fam] = {
            "lengths": lengths,
            "classification": classification,
        }

    # --------------------------------------------------------
    # Historical RFD3 outputs
    # --------------------------------------------------------

    historical_files = []

    for root in (
        slyb / "out",
        slyb / "anansi",
    ):
        if not root.exists():
            continue

        historical_files.extend(
            root.rglob("*.cif.gz")
        )

    # Avoid pathological unlimited scans.
    historical_files = sorted(
        set(historical_files)
    )[:5000]

    print(
        "historical RFD3 structures:",
        len(historical_files),
    )

    historical_chains = (
        defaultdict(list)
    )

    historical_errors = []

    for i, path in enumerate(
        historical_files,
        1,
    ):
        if i % 100 == 0:
            print(
                f"historical scan {i}/"
                f"{len(historical_files)}"
            )

        try:
            data = structure_data(path)

            for c in data:
                historical_chains[
                    c["n_res"]
                ].append(
                    {
                        "path": path,
                        "chain": c["chain"],
                        "seq": c["seq"],
                        "coords": c[
                            "coords"
                        ],
                    }
                )

        except Exception as e:
            historical_errors.append(
                {
                    "path": str(path),
                    "error": repr(e),
                }
            )

    # --------------------------------------------------------
    # Match current parent chains back to historical RFD3
    # --------------------------------------------------------

    hist_matches = []

    for fam, chains in parent_data.items():
        for pc in chains:
            possible = (
                historical_chains.get(
                    pc["n_res"],
                    [],
                )
            )

            scored = []

            for h in possible:
                rmsd = kabsch_rmsd(
                    pc["coords"],
                    h["coords"],
                )

                if rmsd is None:
                    continue

                scored.append(
                    (
                        rmsd,
                        seq_identity(
                            pc["seq"],
                            h["seq"],
                        ),
                        h,
                    )
                )

            scored.sort(
                key=lambda x: x[0]
            )

            for rank, (
                rmsd,
                sid,
                h,
            ) in enumerate(
                scored[:20],
                1,
            ):
                sidecar = sidecar_for(
                    h["path"]
                )

                meta = (
                    load_json(sidecar)
                    if sidecar
                    and sidecar.exists()
                    else None
                )

                target_input = (
                    nested_value(
                        meta,
                        "input",
                    )
                    if meta
                    else None
                )

                specification = (
                    meta.get(
                        "specification"
                    )
                    if isinstance(
                        meta,
                        dict,
                    )
                    else None
                )

                spec_text = (
                    json.dumps(
                        specification,
                        separators=(
                            ",",
                            ":",
                        ),
                    )[:1000]
                    if specification
                    is not None
                    else ""
                )

                hist_matches.append(
                    {
                        "family": fam,
                        "parent_chain": (
                            pc["chain"]
                        ),
                        "rank": rank,
                        "historical_path": (
                            str(
                                h["path"]
                            )
                        ),
                        "historical_chain": (
                            h["chain"]
                        ),
                        "CA_RMSD_A": round(
                            rmsd,
                            6,
                        ),
                        "sequence_identity": (
                            round(
                                sid,
                                6,
                            )
                        ),
                        "sidecar": (
                            str(sidecar)
                            if sidecar
                            else ""
                        ),
                        "target_input": (
                            target_input
                            if isinstance(
                                target_input,
                                str,
                            )
                            else ""
                        ),
                        "specification": (
                            spec_text
                        ),
                    }
                )

    # --------------------------------------------------------
    # Cofold sample inventory + ranking
    # --------------------------------------------------------

    cofold_rows = []

    cofold_root = (
        slyb
        / "cofolding"
        / "out"
    )

    summaries = sorted(
        cofold_root.rglob(
            "*_summary_confidences.json"
        )
    )

    for summary_path in summaries:
        cand = candidate_name(
            summary_path
        )

        if not cand:
            continue

        meta = load_json(
            summary_path
        ) or {}

        model_path = Path(
            str(summary_path).replace(
                "_summary_confidences.json",
                "_model.cif",
            )
        )

        lengths = []
        native_ids = []
        chains = []

        if model_path.exists():
            try:
                data = structure_data(
                    model_path
                )

                for c in data:
                    chains.append(
                        c["chain"]
                    )
                    lengths.append(
                        c["n_res"]
                    )

                    native_ids.append(
                        round(
                            seq_identity(
                                c["seq"],
                                native_a[
                                    "seq"
                                ],
                            ),
                            4,
                        )
                    )

            except Exception:
                pass

        binder_guess = ""

        if cand in candidate_data:
            monomer_len = (
                candidate_data[
                    cand
                ][0]["n_res"]
            )

            matches = [
                chains[i]
                for i, n in enumerate(
                    lengths
                )
                if n == monomer_len
            ]

            if len(matches) == 1:
                binder_guess = (
                    matches[0]
                )

        cofold_rows.append(
            {
                "candidate": cand,
                "summary": str(
                    summary_path
                ),
                "model": str(
                    model_path
                ),
                "chain_count": len(
                    lengths
                ),
                "chains": ",".join(
                    chains
                ),
                "chain_lengths": ",".join(
                    map(str, lengths)
                ),
                "native_seq_identity": (
                    ",".join(
                        map(
                            str,
                            native_ids,
                        )
                    )
                ),
                "binder_chain_guess": (
                    binder_guess
                ),
                "ranking_score": (
                    meta.get(
                        "ranking_score",
                        "",
                    )
                ),
                "iptm": meta.get(
                    "iptm",
                    "",
                ),
                "ptm": meta.get(
                    "ptm",
                    "",
                ),
                "overall_plddt": (
                    meta.get(
                        "overall_plddt",
                        "",
                    )
                ),
                "overall_pae": (
                    meta.get(
                        "overall_pae",
                        "",
                    )
                ),
                "has_clash": (
                    meta.get(
                        "has_clash",
                        "",
                    )
                ),
            }
        )

    # --------------------------------------------------------
    # Existing pair checker on 2-chain structures
    # --------------------------------------------------------

    pair_rows = []

    checker = (
        repo
        / "scripts"
        / "slyb"
        / "check_binder_symmetry.py"
    )

    pair_inputs = list(
        parents.values()
    )

    pair_inputs.extend(
        sorted(
            (
                slyb
                / "sequence_design"
                / "out"
            ).glob(
                "lmpnn_*/*.cif"
            )
        )
    )

    for path in pair_inputs:
        try:
            proc = subprocess.run(
                [
                    sys.executable,
                    str(checker),
                    str(path),
                    "--pair",
                ],
                capture_output=True,
                text=True,
                cwd=slyb,
            )

            obj = json.loads(
                proc.stdout
            )

            first = (
                obj[0]
                if isinstance(
                    obj,
                    list,
                )
                and obj
                else {}
            )

            pair_rows.append(
                {
                    "path": str(path),
                    "exit_code": (
                        proc.returncode
                    ),
                    "PASS": first.get(
                        "PASS",
                        "",
                    ),
                    "binder_binder_min_A": (
                        first.get(
                            "binder_binder_min_A",
                            "",
                        )
                    ),
                    "binder_binder_clash": (
                        first.get(
                            "binder_binder_clash",
                            "",
                        )
                    ),
                    "raw": json.dumps(
                        first,
                        separators=(
                            ",",
                            ":",
                        ),
                    ),
                }
            )

        except Exception as e:
            pair_rows.append(
                {
                    "path": str(path),
                    "exit_code": "",
                    "PASS": "",
                    "binder_binder_min_A": "",
                    "binder_binder_clash": "",
                    "raw": repr(e),
                }
            )

    # --------------------------------------------------------
    # Write tables
    # --------------------------------------------------------

    write_tsv(
        out / "structures.tsv",
        structure_rows,
        [
            "group",
            "name",
            "path",
            "chain",
            "n_res",
            "seq_sha1",
            "native_seq_identity",
        ],
    )

    write_tsv(
        out
        / "parent_monomer_matches.tsv",
        monomer_matches,
        [
            "candidate",
            "family",
            "parent_chain",
            "parent_n_res",
            "monomer_n_res",
            "sequence_identity",
            "CA_RMSD_A",
        ],
    )

    write_tsv(
        out
        / "historical_rfd3_matches.tsv",
        hist_matches,
        [
            "family",
            "parent_chain",
            "rank",
            "historical_path",
            "historical_chain",
            "CA_RMSD_A",
            "sequence_identity",
            "sidecar",
            "target_input",
            "specification",
        ],
    )

    write_tsv(
        out / "cofold_samples.tsv",
        cofold_rows,
        [
            "candidate",
            "summary",
            "model",
            "chain_count",
            "chains",
            "chain_lengths",
            "native_seq_identity",
            "binder_chain_guess",
            "ranking_score",
            "iptm",
            "ptm",
            "overall_plddt",
            "overall_pae",
            "has_clash",
        ],
    )

    write_tsv(
        out / "pair_check.tsv",
        pair_rows,
        [
            "path",
            "exit_code",
            "PASS",
            "binder_binder_min_A",
            "binder_binder_clash",
            "raw",
        ],
    )

    write_tsv(
        out / "historical_errors.tsv",
        historical_errors,
        [
            "path",
            "error",
        ],
    )

    # --------------------------------------------------------
    # Best cofold per candidate
    # --------------------------------------------------------

    best_cofold = {}

    for cand in candidates:
        rows = [
            x
            for x in cofold_rows
            if x["candidate"] == cand
        ]

        def score(row):
            clash = str(
                row["has_clash"]
            ).lower() == "true"

            try:
                ranking = float(
                    row[
                        "ranking_score"
                    ]
                )
            except Exception:
                ranking = -999

            try:
                iptm = float(
                    row["iptm"]
                )
            except Exception:
                iptm = -999

            return (
                not clash,
                ranking,
                iptm,
            )

        if rows:
            best_cofold[cand] = max(
                rows,
                key=score,
            )

    # --------------------------------------------------------
    # JSON summary
    # --------------------------------------------------------

    top_hist = defaultdict(list)

    for row in hist_matches:
        if row["rank"] <= 5:
            top_hist[
                f"{row['family']}:{row['parent_chain']}"
            ].append(row)

    summary = {
        "native_chain_A_residues": (
            native_a["n_res"]
        ),
        "parent_classification": (
            parent_class
        ),
        "candidate_family": family,
        "best_cofold": best_cofold,
        "top_historical_matches": dict(
            top_hist
        ),
        "historical_structure_count": (
            len(historical_files)
        ),
        "cofold_sample_count": (
            len(cofold_rows)
        ),
        "pair_check_count": (
            len(pair_rows)
        ),
    }

    (
        out
        / "lineage_summary.json"
    ).write_text(
        json.dumps(
            summary,
            indent=2,
        )
        + "\n"
    )

    # --------------------------------------------------------
    # Human-readable summary
    # --------------------------------------------------------

    lines = []

    lines.append(
        "# SlyB lineage recovery"
    )

    lines.append("")

    lines.append(
        f"Native SlyB chain A: "
        f"{native_a['n_res']} residues"
    )

    lines.append(
        f"Historical RFD3 structures scanned: "
        f"{len(historical_files)}"
    )

    lines.append(
        f"Cofold samples indexed: "
        f"{len(cofold_rows)}"
    )

    lines.append("")

    lines.append(
        "## Parent representations"
    )

    lines.append("")

    for fam, info in parent_class.items():
        lines.append(
            f"- {fam}: chains "
            f"{info['lengths']} -> "
            f"**{info['classification']}**"
        )

    lines.append("")

    lines.append(
        "## RF3 monomer -> parent matches"
    )

    lines.append("")

    for cand in candidates:
        rows = [
            x
            for x in monomer_matches
            if x["candidate"] == cand
        ]

        rows.sort(
            key=lambda x: (
                999999
                if x["CA_RMSD_A"] == ""
                else float(
                    x["CA_RMSD_A"]
                )
            )
        )

        if rows:
            x = rows[0]

            lines.append(
                f"- {cand} -> "
                f"{x['family']} chain "
                f"{x['parent_chain']}: "
                f"CA RMSD={x['CA_RMSD_A']} Å, "
                f"seq identity="
                f"{x['sequence_identity']}"
            )

    lines.append("")

    lines.append(
        "## Best cofold samples"
    )

    lines.append("")

    for cand, x in best_cofold.items():
        lines.append(
            f"- {cand}: "
            f"ranking={x['ranking_score']}, "
            f"iPTM={x['iptm']}, "
            f"chains={x['chain_lengths']}, "
            f"binder={x['binder_chain_guess']}, "
            f"clash={x['has_clash']}"
        )

    lines.append("")

    lines.append(
        "## Next decision"
    )

    lines.append("")

    lines.append(
        "Use historical_rfd3_matches.tsv "
        "to restore the original target/input "
        "specification before authorizing new RFD3."
    )

    (
        out / "SUMMARY.md"
    ).write_text(
        "\n".join(lines)
        + "\n"
    )

    print()
    print(
        "\n".join(lines)
    )


if __name__ == "__main__":
    main()
