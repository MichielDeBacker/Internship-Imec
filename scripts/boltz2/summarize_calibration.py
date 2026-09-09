from pathlib import Path
import csv
import json
import statistics

ROOT = Path.home() / "my_project"

roles = {}
with (ROOT / "configs/boltz2/binder_set.tsv").open() as f:
    for row in csv.DictReader(f, delimiter="\t"):
        roles[row["name"]] = row["role"]

def pair(d, i, j):
    p = d.get("pair_chains_iptm", {})
    try:
        return float(p[str(i)][str(j)])
    except (KeyError, TypeError, ValueError):
        return None

def mean_present(*values):
    vals = [v for v in values if v is not None]
    return statistics.mean(vals) if vals else None

records = []

# Existing native positive control.
native = list(
    (ROOT / "results/boltz2_calibration/native_AB").rglob(
        "confidence_native_AB_model_0.json"
    )
)

if native:
    d = json.loads(native[0].read_text())
    records.append({
        "name": "native_AB",
        "role": "positive_control",
        "iptm": d.get("iptm"),
        "protein_iptm": d.get("protein_iptm"),
        "complex_plddt": d.get("complex_plddt"),
        "AB": mean_present(pair(d, 0, 1), pair(d, 1, 0)),
        "AC": None,
        "BC": None,
        "binder_interface": None,
    })

# Candidate + negative-control complexes.
base = ROOT / "results/boltz2_calibration/all"

for conf in sorted(base.rglob("confidence_*_model_0.json")):
    name = conf.name.removeprefix("confidence_").removesuffix("_model_0.json")
    d = json.loads(conf.read_text())

    ab = mean_present(pair(d, 0, 1), pair(d, 1, 0))
    ac = mean_present(pair(d, 0, 2), pair(d, 2, 0))
    bc = mean_present(pair(d, 1, 2), pair(d, 2, 1))

    records.append({
        "name": name,
        "role": roles.get(name, "unknown"),
        "iptm": d.get("iptm"),
        "protein_iptm": d.get("protein_iptm"),
        "complex_plddt": d.get("complex_plddt"),
        "AB": ab,
        "AC": ac,
        "BC": bc,
        "binder_interface": mean_present(ac, bc),
    })

def fmt(x):
    return "" if x is None else f"{float(x):.4f}"

tsv = ROOT / "reports/boltz2/calibration_summary.tsv"

with tsv.open("w") as f:
    f.write(
        "name\trole\tiptm\tprotein_iptm\tcomplex_plddt\t"
        "AB_pair_iptm\tAC_pair_iptm\tBC_pair_iptm\tbinder_interface_mean\n"
    )
    for r in records:
        f.write(
            "\t".join([
                r["name"],
                r["role"],
                fmt(r["iptm"]),
                fmt(r["protein_iptm"]),
                fmt(r["complex_plddt"]),
                fmt(r["AB"]),
                fmt(r["AC"]),
                fmt(r["BC"]),
                fmt(r["binder_interface"]),
            ]) + "\n"
        )

print(tsv.read_text())
