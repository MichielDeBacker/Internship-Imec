#!/bin/bash
#SBATCH --job-name=slyb_rfd3
#SBATCH --account=bvo00014
#SBATCH --partition=ampere_gpu
#SBATCH --gpus=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# RFdiffusion3 groove-targeted staple design on hydra.
#   partition ampere_gpu : A100 (2/node, 256 GB, 32 cores) -- measured, not assumed
#   defaults per node    : 8 cores + 64 GB is a quarter node per GPU
#   MaxTime on partition : 5-00:00:00, so 8h is well inside
#
# Submit:  sbatch run_rfd3_hydra.sh                  (all arms)
#          sbatch run_rfd3_hydra.sh groove_80_duo    (one arm)
set -euo pipefail

ROOT="${SLYB_ROOT:-/scratch/brussel/vo/000/bvo00014/vsc39230/slyb}"
ARM="${1:-all}"

# ---- preflight 0: the Slurm association must actually be the new VO --------
# The Unix group switched to bvo00014 before the Slurm association did. A job
# submitted while sacctmgr still reports bvo00004 would bill the wrong VO, so
# refuse rather than silently charge the old allocation.
WANT_ACCT="bvo00014"
HAVE_ACCTS=$(sacctmgr -nP show assoc user="$USER" format=Account 2>/dev/null | sort -u | tr '\n' ' ')
echo "=== slurm associations: ${HAVE_ACCTS:-none} ==="
case " $HAVE_ACCTS " in
    *" $WANT_ACCT "*) echo "account $WANT_ACCT present -- ok" ;;
    *) echo "FATAL: sacctmgr does not list $WANT_ACCT for $USER (has: ${HAVE_ACCTS:-none})."
       echo "       The VO association has not propagated yet. Not submitting."
       exit 1 ;;
esac

# ---- preflight: refuse to burn GPU hours on a broken environment ----------
if [ ! -f "$ROOT/env.sh" ]; then
    echo "FATAL: $ROOT/env.sh missing -- run install_hydra.sh on the login node first"
    exit 1
fi
# shellcheck disable=SC1091
source "$ROOT/env.sh"

echo "=== preflight ==="
if ! python "$ROOT/verify_hydra.py"; then
    echo "FATAL: verification failed -- not launching. Fix the FAIL rows above."
    exit 1
fi

INPUTS="$ROOT/inputs/rfd3_groove_inputs.json"
[ -f "$INPUTS" ] || { echo "FATAL: inputs missing: $INPUTS"; exit 1; }

for f in target_duo_soluble.pdb target_trio_soluble.pdb; do
    [ -f "$ROOT/inputs/$f" ] || { echo "FATAL: target missing: $ROOT/inputs/$f"; exit 1; }
done

# input paths inside the JSON must resolve on THIS filesystem
python - "$INPUTS" <<'PYCHK'
import json, os, sys
spec = json.load(open(sys.argv[1]))
missing = [(k, v["input"]) for k, v in spec.items() if not os.path.isfile(v["input"])]
if missing:
    for k, p in missing:
        print("FATAL: arm %s points at a non-existent input: %s" % (k, p))
    sys.exit(1)
print("all %d arms have resolvable inputs" % len(spec))
PYCHK

echo "=== dry run (prevalidate_inputs) ==="
rfd3 design inputs="$INPUTS" out_dir="$ROOT/out/_prevalidate" \
    prevalidate_inputs=True
echo "spec validation passed"

# ---- production ----------------------------------------------------------
OUT="$ROOT/out/$ARM"
mkdir -p "$OUT"

echo "=== GPU at launch ==="
nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv,noheader

# A100 80GB: batch 4 is comfortable; low_memory_mode deliberately OFF because
# it held utilisation near 14% in the reference benchmark.
echo "=== rfd3 design: arm=$ARM ==="
if [ "$ARM" = "all" ]; then
    rfd3 design inputs="$INPUTS" out_dir="$OUT" \
        n_batches=8 diffusion_batch_size=4
else
    python - "$INPUTS" "$ARM" "$ROOT/inputs/_arm_${ARM}.json" <<'PYCHK'
import json, sys
allspec = json.load(open(sys.argv[1]))
arm = sys.argv[2]
if arm not in allspec:
    sys.exit("FATAL: arm %r not in spec. Available: %s" % (arm, ", ".join(allspec)))
json.dump({arm: allspec[arm]}, open(sys.argv[3], "w"), indent=2)
print("wrote single-arm spec for", arm)
PYCHK
    rfd3 design inputs="$ROOT/inputs/_arm_${ARM}.json" out_dir="$OUT" \
        n_batches=8 diffusion_batch_size=4
fi

echo "=== designs produced ==="
ls -1 "$OUT"/*.cif.gz 2>/dev/null | wc -l

echo "=== symmetry / steric triage ==="
python "$ROOT/check_binder_symmetry.py" "$OUT"/*.cif.gz \
    --ref "$ROOT/inputs/target_ring_C11.pdb" > "$OUT/symmetry_report.json" || \
    echo "WARN: triage step failed -- designs are kept, re-run triage manually"

echo "=== done: $(date -Is) ==="
