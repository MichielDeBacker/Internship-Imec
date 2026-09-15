#!/usr/bin/env python3

import argparse
import csv
import hashlib
import importlib
import importlib.metadata
import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "pipeline" / "project.json"
WORK = ROOT / "work"
REPORTS = ROOT / "reports" / "slyb"


def load_config():
    with CONFIG_PATH.open() as f:
        return json.load(f)


CFG = load_config()


def slyb_root():
    return Path(
        os.environ.get(
            "SLYB_ROOT",
            CFG["default_slyb_root"]
        )
    ).resolve()


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def rel(path, root):
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def git_branch():
    try:
        return subprocess.check_output(
            ["git", "branch", "--show-current"],
            cwd=ROOT,
            text=True,
        ).strip()
    except Exception:
        return "UNKNOWN"


def candidate_dict():
    out = {}
    out.update(CFG["primary_candidates"])
    out.update(CFG["reserve_and_controls"])
    return out


def add_record(records, root, path, role):
    if not path.exists() or not path.is_file():
        return

    records.append(
        {
            "role": role,
            "path": rel(path, root),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    )


def cmd_status(_args=None):
    root = slyb_root()

    print("project       :", CFG["project"])
    print("repo          :", ROOT)
    print("git branch    :", git_branch())
    print("SLYB_ROOT     :", root)
    print("work          :", WORK)

    print("\nTargets:")
    for name, value in CFG["targets"].items():
        p = root / value
        print(f"  {name:10s} {'OK' if p.exists() else 'MISSING'}  {p}")

    print("\nBackbone families:")
    for name, value in CFG["designs"].items():
        p = root / value
        print(f"  {name:10s} {'OK' if p.exists() else 'MISSING'}  {p}")

    print("\nCandidates:")
    for name, value in candidate_dict().items():
        p = root / value
        print(f"  {name:10s} {'OK' if p.exists() else 'MISSING'}  {p}")

    inv = WORK / "recovered_inventory.json"
    geo = WORK / "geometry" / "summary.tsv"
    plan = WORK / "rfd3" / "plan.json"

    print("\nPipeline artifacts:")
    print("  inventory :", "READY" if inv.exists() else "NOT RUN")
    print("  geometry  :", "READY" if geo.exists() else "NOT RUN")
    print("  rfd3 plan :", "READY" if plan.exists() else "NOT RUN")

    print("\nRFD3 runtime:")
    print("  rfd3      :", shutil.which("rfd3") or "NOT ON PATH")
    print("  foundry   :", shutil.which("foundry") or "NOT ON PATH")

    return 0


def cmd_doctor(_args=None):
    root = slyb_root()

    print("=== Python ===")
    print(sys.executable)
    print(sys.version.replace("\n", " "))

    ok = True

    print("\n=== Python dependencies ===")
    for module in ("numpy", "scipy", "biotite"):
        try:
            m = importlib.import_module(module)
            version = getattr(m, "__version__", "unknown")
            print(f"OK      {module:10s} {version}")
        except Exception as e:
            print(f"MISSING {module:10s} {e}")
            ok = False

    print("\n=== Structural inputs ===")

    required = []

    for name, value in CFG["targets"].items():
        required.append((f"target:{name}", root / value))

    for name, value in CFG["designs"].items():
        required.append((f"design:{name}", root / value))

    required.append(
        (
            "geometry-script",
            ROOT / "scripts" / "slyb" / "check_binder_symmetry.py",
        )
    )

    for role, path in required:
        exists = path.is_file()
        print(
            f"{'OK' if exists else 'MISSING':7s} "
            f"{role:20s} {path}"
        )
        ok = ok and exists

    return 0 if ok else 2


def cmd_inventory(_args=None):
    root = slyb_root()
    WORK.mkdir(parents=True, exist_ok=True)

    records = []

    for name, value in CFG["targets"].items():
        add_record(records, root, root / value, f"target:{name}")

    for name, value in CFG["designs"].items():
        add_record(records, root, root / value, f"backbone:{name}")

    for name, value in candidate_dict().items():
        model = root / value
        add_record(records, root, model, f"candidate:{name}")

        parent = model.parent

        for suffix in (
            "_confidences.json",
            "_summary_confidences.json",
            "_ranking_scores.csv",
        ):
            add_record(
                records,
                root,
                parent / f"{name}{suffix}",
                f"candidate-metadata:{name}",
            )

    cofold_root = root / "cofolding" / "out"

    if cofold_root.exists():
        for name in candidate_dict():
            for p in sorted(cofold_root.rglob(f"*{name}*")):
                if not p.is_file():
                    continue

                low = p.name.lower()

                if (
                    low.endswith(".cif")
                    or low.endswith(".cif.gz")
                    or low.endswith(".json")
                    or low.endswith(".csv")
                    or low.endswith(".npz")
                ):
                    add_record(
                        records,
                        root,
                        p,
                        f"cofold:{name}",
                    )

    out = WORK / "recovered_inventory.json"

    payload = {
        "generated_at": now(),
        "project": CFG["project"],
        "slyb_root": str(root),
        "record_count": len(records),
        "records": records,
    }

    with out.open("w") as f:
        json.dump(payload, f, indent=2)

    print(f"wrote {out}")
    print(f"records: {len(records)}")

    roles = Counter(r["role"].split(":")[0] for r in records)

    for key, value in sorted(roles.items()):
        print(f"  {key:20s} {value}")

    return 0


def parse_pass(text):
    patterns = [
        r"""["']?PASS["']?\s*[:=]\s*(True|False|true|false)""",
        r"""\bPASS\s+(True|False|true|false)\b""",
    ]

    for pattern in patterns:
        m = re.search(pattern, text)
        if m:
            return "PASS" if m.group(1).lower() == "true" else "FAIL"

    return "UNKNOWN"


def cmd_geometry(_args=None):
    root = slyb_root()

    ring = root / CFG["targets"]["ring"]
    script = ROOT / "scripts" / "slyb" / "check_binder_symmetry.py"

    outdir = WORK / "geometry"
    logs = outdir / "logs"

    logs.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    structures = []

    for name, value in CFG["designs"].items():
        p = root / value
        if p.is_file():
            structures.append(
                {
                    "name": name,
                    "class": "original_backbone",
                    "path": p,
                }
            )

    seq_out = root / "sequence_design" / "out"

    if seq_out.exists():
        for p in sorted(seq_out.glob("lmpnn_*/*.cif")):
            structures.append(
                {
                    "name": p.stem,
                    "class": "ligandmpnn",
                    "path": p,
                }
            )

    print(f"structures to screen: {len(structures)}")

    rows = []

    for i, item in enumerate(structures, 1):
        path = item["path"]

        safe = re.sub(
            r"[^A-Za-z0-9_.-]+",
            "_",
            rel(path, root),
        )

        logfile = logs / f"{safe}.txt"

        cmd = [
            sys.executable,
            str(script),
            str(path),
            "--order",
            "11",
            "--ref",
            str(ring),
        ]

        print(
            f"[{i:03d}/{len(structures):03d}] "
            f"{item['class']} {rel(path, root)}"
        )

        proc = subprocess.run(
            cmd,
            cwd=root,
            text=True,
            capture_output=True,
        )

        text = proc.stdout + "\n" + proc.stderr

        logfile.write_text(text)

        status = (
            "ERROR"
            if proc.returncode != 0
            else parse_pass(text)
        )

        rows.append(
            {
                "name": item["name"],
                "class": item["class"],
                "status": status,
                "exit_code": proc.returncode,
                "path": rel(path, root),
                "log": str(logfile),
            }
        )

    summary_tsv = outdir / "summary.tsv"

    with summary_tsv.open("w", newline="") as f:
        fields = [
            "name",
            "class",
            "status",
            "exit_code",
            "path",
            "log",
        ]

        writer = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=fields,
        )

        writer.writeheader()
        writer.writerows(rows)

    summary_json = outdir / "summary.json"

    summary_json.write_text(
        json.dumps(
            {
                "generated_at": now(),
                "rows": rows,
            },
            indent=2,
        )
    )

    shutil.copy2(
        summary_tsv,
        REPORTS / "latest_geometry_summary.tsv",
    )

    counts = Counter(r["status"] for r in rows)

    print("\nGeometry results:")

    for key in ("PASS", "FAIL", "UNKNOWN", "ERROR"):
        print(f"  {key:8s} {counts.get(key, 0)}")

    print(f"\nsummary: {summary_tsv}")

    return 0 if counts.get("ERROR", 0) == 0 else 3


def cmd_rfd3_check(_args=None):
    print("=== RFD3 runtime ===")

    print("rfd3   :", shutil.which("rfd3") or "NOT ON PATH")
    print("foundry:", shutil.which("foundry") or "NOT ON PATH")

    print(
        "FOUNDRY_CHECKPOINT_DIRS:",
        os.environ.get(
            "FOUNDRY_CHECKPOINT_DIRS",
            "NOT SET",
        ),
    )

    print("\nInstalled packages:")

    for package in (
        "rc-foundry",
        "rfd3",
    ):
        try:
            print(
                f"{package}: "
                f"{importlib.metadata.version(package)}"
            )
        except importlib.metadata.PackageNotFoundError:
            print(f"{package}: NOT INSTALLED")

    return 0


def cmd_rfd3_plan(_args=None):
    summary = WORK / "geometry" / "summary.tsv"

    if not summary.exists():
        print(
            "Geometry summary not found. "
            "Run geometry first.",
            file=sys.stderr,
        )
        return 2

    rows = []

    with summary.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows.extend(reader)

    viable = [
        {
            "name": r["name"],
            "class": r["class"],
            "path": r["path"],
        }
        for r in rows
        if r["status"] == "PASS"
    ]

    outdir = WORK / "rfd3"
    outdir.mkdir(parents=True, exist_ok=True)

    payload = {
        "generated_at": now(),
        "authorized": False,
        "reason": (
            "Manual parent selection and target hotspot "
            "review are required before RFD3 execution."
        ),
        "geometry_pass_candidates": viable,
        "maximum_parent_geometries_after_review": 2,
        "hotspots": None,
        "rfd3_command": None,
    }

    out = outdir / "plan.json"

    out.write_text(json.dumps(payload, indent=2))

    REPORTS.mkdir(parents=True, exist_ok=True)

    shutil.copy2(
        out,
        REPORTS / "latest_rfd3_gate.json",
    )

    print(f"wrote {out}")
    print(f"geometry PASS candidates: {len(viable)}")
    print("RFD3 execution authorized: NO")
    print("Reason: hotspot/parent review still required")

    return 0


def cmd_bootstrap(_args=None):
    rc = cmd_doctor()

    if rc != 0:
        print("\nDoctor failed; bootstrap stopped.")
        return rc

    print("\n=== INVENTORY ===")
    rc = cmd_inventory()

    if rc != 0:
        return rc

    print("\n=== GEOMETRY ===")
    rc = cmd_geometry()

    print("\n=== STATUS ===")
    cmd_status()

    return rc


def main():
    parser = argparse.ArgumentParser(
        description=(
            "SlyB reuse-first binder design pipeline "
            "for VUB Hydra."
        )
    )

    sub = parser.add_subparsers(
        dest="command",
        required=True,
    )

    commands = {
        "status": cmd_status,
        "doctor": cmd_doctor,
        "inventory": cmd_inventory,
        "geometry": cmd_geometry,
        "rfd3-check": cmd_rfd3_check,
        "rfd3-plan": cmd_rfd3_plan,
        "bootstrap": cmd_bootstrap,
    }

    for name, func in commands.items():
        p = sub.add_parser(name)
        p.set_defaults(func=func)

    args = parser.parse_args()

    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
