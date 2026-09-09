from pathlib import Path
import csv
import json
import os

ROOT = Path.home() / "my_project"

RESULTS = (
    Path("/emdata/EMprocessing")
    / os.environ["USER"]
    / "slyb_boltz2"
    / "results"
    / "free96_cofolding"
)

OUT = (
    ROOT
    / "reports"
    / "boltz2"
    / "free96_cofolding_summary.tsv"
)

roles = {}

with (
    ROOT
    / "configs"
    / "boltz2"
    / "binder_set.tsv"
).open() as f:
    for row in csv.DictReader(f, delimiter="\t"):
        roles[row["name"]] = row["role"]


def pair(d, a, b):
    p = d.get("pair_chains_iptm", {})
    vals = []

    for x, y in [(a, b), (b, a)]:
        try:
            vals.append(float(p[str(x)][str(y)]))
        except Exception:
            pass

    return sum(vals) / len(vals) if vals else None


rows = []

for path in sorted(
    RESULTS.rglob("confidence_*_model_0.json")
):
    name = (
        path.name
        .removeprefix("confidence_")
        .removesuffix("_model_0.json")
    )

    if name not in roles:
        continue

    d = json.loads(path.read_text())

    ab = pair(d, 0, 1)
    ac = pair(d, 0, 2)
    bc = pair(d, 1, 2)

    if None in (ab, ac, bc):
        raise RuntimeError(
            f"{name}: missing pairwise iPTM"
        )

    rows.append({
        "name": name,
        "role": roles[name],
        "iptm": float(d["iptm"]),
        "AB": ab,
        "AC": ac,
        "BC": bc,
        "mean": (ac + bc) / 2,
        "minimum": min(ac, bc),
        "delta": abs(ac - bc),
    })

rows.sort(
    key=lambda r: -r["minimum"]
)

OUT.parent.mkdir(
    parents=True,
    exist_ok=True,
)

with OUT.open(
    "w",
    newline="",
) as f:
    w = csv.writer(
        f,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writerow([
        "name",
        "role",
        "iptm",
        "AB_pair_iptm",
        "AC_pair_iptm",
        "BC_pair_iptm",
        "binder_mean",
        "binder_min",
        "binder_delta",
    ])

    for r in rows:
        w.writerow([
            r["name"],
            r["role"],
            f'{r["iptm"]:.4f}',
            f'{r["AB"]:.4f}',
            f'{r["AC"]:.4f}',
            f'{r["BC"]:.4f}',
            f'{r["mean"]:.4f}',
            f'{r["minimum"]:.4f}',
            f'{r["delta"]:.4f}',
        ])

print(OUT.read_text())
