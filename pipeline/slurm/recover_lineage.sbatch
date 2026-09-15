#!/bin/bash
#SBATCH --job-name=slyb_lineage
#SBATCH --partition=zen4
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:55:00
#SBATCH --output=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/lineage_%j.out
#SBATCH --error=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/lineage_%j.out

set -uo pipefail

BASE=/scratch/brussel/vo/000/bvo00014/vsc39230
SLYB="$BASE/slyb"
REPO="$BASE/repos/Internship-Imec"
ENV="$BASE/venvs/slyb-geometry"

source "$ENV/bin/activate"

cd "$REPO"

OUT="$REPO/work/lineage_recovery/${SLURM_JOB_ID}"

mkdir -p \
    "$OUT" \
    "$REPO/reports/slyb/lineage_latest"

echo "============================================================"
echo "SLYB LINEAGE RECOVERY"
echo "job=$SLURM_JOB_ID"
echo "host=$(hostname)"
echo "start=$(date)"
echo "out=$OUT"
echo "============================================================"

echo
echo "=== PYTHON ENV ==="

python --version

python - <<'PY'
import numpy
import scipy
import biotite

print("numpy", numpy.__version__)
print("scipy", scipy.__version__)
print("biotite", biotite.__version__)
PY

echo
echo "=== PIPELINE STATUS BEFORE ==="

export SLYB_ROOT="$SLYB"

python pipeline/run.py status \
    > "$OUT/pipeline_status_before.txt" \
    2>&1 || true

cat "$OUT/pipeline_status_before.txt"

echo
echo "=== FULL LINEAGE RECOVERY ==="

python pipeline/lib/recover_lineage.py \
    --slyb "$SLYB" \
    --repo "$REPO" \
    --out "$OUT"

echo
echo "=== HISTORICAL TWO-BINDER / PENTA SCREEN ==="

PENTA_ROOT="$SLYB/out/rfd3_T200_penta_single_screen"

if [ -d "$PENTA_ROOT" ]; then

    echo "Found historical penta screen:"
    echo "$PENTA_ROOT"

    python \
        "$REPO/scripts/slyb/check_penta_2b_backbone.py" \
        --root "$PENTA_ROOT" \
        --json-out "$OUT/penta_2binder.json" \
        --csv-out "$OUT/penta_2binder.csv" \
        > "$OUT/penta_checker.log" \
        2>&1 || true

    cat "$OUT/penta_checker.log"

else

    echo "Historical penta directory not found."

fi

echo
echo "=== EXISTING RFD3 SYMMETRY CHECKER HELP ==="

python \
    "$REPO/scripts/slyb/check_rfd3_symmetry.py" \
    --help \
    > "$OUT/check_rfd3_symmetry_help.txt" \
    2>&1 || true

cat "$OUT/check_rfd3_symmetry_help.txt"

echo
echo "=== COPY COMPACT REPORTS INTO REPO ==="

for f in \
    SUMMARY.md \
    lineage_summary.json \
    structures.tsv \
    parent_monomer_matches.tsv \
    historical_rfd3_matches.tsv \
    cofold_samples.tsv \
    pair_check.tsv \
    historical_errors.tsv \
    pipeline_status_before.txt \
    penta_2binder.json \
    penta_2binder.csv \
    penta_checker.log \
    check_rfd3_symmetry_help.txt
do
    if [ -f "$OUT/$f" ]; then
        cp -f \
            "$OUT/$f" \
            "$REPO/reports/slyb/lineage_latest/$f"
    fi
done

echo
echo "=== RESULT COUNTS ==="

for f in \
    structures.tsv \
    parent_monomer_matches.tsv \
    historical_rfd3_matches.tsv \
    cofold_samples.tsv \
    pair_check.tsv
do
    if [ -f "$OUT/$f" ]; then
        echo "$f: $(($(wc -l < "$OUT/$f") - 1)) rows"
    fi
done

echo
echo "=== TOP HISTORICAL MATCHES ==="

if [ -f "$OUT/historical_rfd3_matches.tsv" ]; then
    {
        head -1 \
            "$OUT/historical_rfd3_matches.tsv"

        awk -F '\t' \
            'NR>1 && $3<=3' \
            "$OUT/historical_rfd3_matches.tsv"
    } \
    | column -t -s $'\t' \
    | head -100
fi

echo
echo "=== BEST COFOLD SUMMARY ==="

cat "$OUT/SUMMARY.md"

echo
echo "=== LOCAL GIT STATE ==="

git status --short

echo
echo "============================================================"
echo "LINEAGE RECOVERY FINISHED"
echo "end=$(date)"
echo "OUT=$OUT"
echo "============================================================"

touch "$OUT/DONE"
