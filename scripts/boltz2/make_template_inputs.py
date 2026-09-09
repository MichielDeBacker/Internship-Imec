from pathlib import Path
import csv
import os

ROOT = Path.home() / "my_project"

TEMPLATE = (
    Path("/emdata/EMprocessing")
    / os.environ["USER"]
    / "slyb_boltz2"
    / "data"
    / "7ojg.cif"
)

SLYB = (
    ROOT
    / "configs/boltz2/7ojg_sequence.txt"
).read_text().strip()

assert len(SLYB) == 155, len(SLYB)

table = ROOT / "configs/boltz2/binder_set.tsv"
outdir = ROOT / "configs/boltz2/template_inputs"

outdir.mkdir(
    parents=True,
    exist_ok=True,
)

with table.open() as f:
    rows = list(
        csv.DictReader(
            f,
            delimiter="\t",
        )
    )

designs = [
    (
        r["name"],
        r["role"],
        r["sequence"].strip(),
    )
    for r in rows
]

designs.append(
    (
        "polyA7",
        "sanity_negative",
        "AAAAAAA",
    )
)

for name, role, binder in designs:

    text = f"""version: 1

sequences:
  - protein:
      id: [A, B]
      sequence: {SLYB}
      msa: empty

  - protein:
      id: C
      sequence: {binder}
      msa: empty

templates:
  - cif: {TEMPLATE}
    chain_id: [A, B]
    template_id: [A, B]
    force: true
    threshold: 1.0
"""

    path = outdir / f"{name}.yaml"
    path.write_text(text)

    print(
        f"{name:8s} "
        f"{role:16s} "
        f"A/B={len(SLYB)} aa "
        f"C={len(binder)} aa"
    )
