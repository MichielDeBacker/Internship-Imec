from pathlib import Path
import csv

ROOT = Path.home() / "my_project"

SLYB = (
    "CVNNDTLSGDVYTASEAKQVQNVSYGTIVNVRPVQIQGGD"
    "QGVQSAMNKTQGVELEIRKDDGNTIMVVQKQGNTRFSPGQRVVLASNGSQVTVSPR"
)

assert len(SLYB) == 96, len(SLYB)

table = ROOT / "configs/boltz2/binder_set.tsv"
outdir = ROOT / "configs/boltz2/free96_inputs"
outdir.mkdir(parents=True, exist_ok=True)

with table.open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

expected = {"d1_s3", "d2_s0", "d2_s1", "d6_s3", "d6_s7"}

found = {r["name"] for r in rows}

missing = expected - found
if missing:
    raise RuntimeError(f"Missing binders: {sorted(missing)}")

for row in rows:
    name = row["name"]

    if name not in expected:
        continue

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

    if "templates:" in text:
        raise RuntimeError("Template unexpectedly present")

    path = outdir / f"{name}.yaml"
    path.write_text(text)

    print(
        f"{name:6s} "
        f"{row['role']:10s} "
        f"A=96 B=96 C={len(binder)}"
    )
