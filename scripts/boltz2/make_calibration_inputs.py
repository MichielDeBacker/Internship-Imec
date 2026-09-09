from pathlib import Path
import csv

ROOT = Path.home() / "my_project"

# Exact 96-aa soluble SlyB construct used for the calibration.
SLYB = (
    "CVNNDTLSGDVYTASEAKQVQNVSYGTIVNVRPVQIQGGD"
    "QGVQSAMNKTQGVELEIRKDDGNTIMVVQKQGNTRFSPGQRVVLASNGSQVTVSPR"
)

assert len(SLYB) == 96

table = ROOT / "configs/boltz2/binder_set.tsv"
outdir = ROOT / "configs/boltz2/calibration_inputs"
outdir.mkdir(parents=True, exist_ok=True)

with table.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

for row in rows:
    name = row["name"]
    binder = row["sequence"].strip()

    text = f"""version: 1

sequences:
  - protein:
      id: A
      sequence: {SLYB}
      msa: empty

  - protein:
      id: B
      sequence: {SLYB}
      msa: empty

  - protein:
      id: C
      sequence: {binder}
      msa: empty
"""

    path = outdir / f"{name}.yaml"
    path.write_text(text)
    print(f"{name:8s} {row['role']:10s} binder_len={len(binder):2d} -> {path}")
