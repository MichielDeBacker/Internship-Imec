#!/bin/bash
#SBATCH --job-name=slyb_vertical
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/vertical_rfd3_%j.out
#SBATCH --error=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/vertical_rfd3_%j.out

set -Eeuo pipefail

BASE=/scratch/brussel/vo/000/bvo00014/vsc39230
SLYB="$BASE/slyb"
REPO="$BASE/repos/Internship-Imec"
GEOMENV="$BASE/venvs/slyb-geometry"
RFDENV="$BASE/venvs/rfd3-foundry"
CKPT="$BASE/models/foundry_checkpoints"
RUN="$SLYB/out/rfd3_vertical_${SLURM_JOB_ID}"
TARGET="$RUN/target_AB_vertical.pdb"
SPEC="$RUN/vertical_spec.json"
GEO="$RUN/vertical_geometry.json"
PRE="$RUN/preflight"
PROD="$RUN/production"
SCREEN="$RUN/screen"

mkdir -p "$RUN" "$PRE" "$PROD" "$SCREEN" "$CKPT" "$SLYB/logs"
cd "$REPO"

echo "============================================================"
echo "SLYB VERTICAL-GROOVE RFD3"
echo "JOB=$SLURM_JOB_ID"
echo "HOST=$(hostname)"
echo "START=$(date)"
echo "RUN=$RUN"
echo "============================================================"

nvidia-smi

source "$GEOMENV/bin/activate"
python pipeline/lib/build_vertical_groove_spec.py \
  --ring "$SLYB/inputs/target_ring_C11.pdb" \
  --penta "$SLYB/inputs/target_penta_soluble.pdb" \
  --parent "$SLYB/sequence_design/inputs/design2.cif" \
  --target-out "$TARGET" \
  --spec-out "$SPEC" \
  --geometry-out "$GEO"

cat "$GEO"

deactivate || true
if [ ! -x "$RFDENV/bin/rfd3" ]; then
  python3 -m venv "$RFDENV"
  source "$RFDENV/bin/activate"
  python -m pip install --upgrade pip setuptools wheel
  python -m pip install "rc-foundry[rfd3]"
else
  source "$RFDENV/bin/activate"
fi

export FOUNDRY_CHECKPOINT_DIRS="$CKPT"
if ! find "$CKPT" -type f -name '*.ckpt' -print -quit | grep -q .; then
  foundry install rfd3 --checkpoint-dir "$CKPT"
fi

python - <<'PY'
import torch
print('CUDA:', torch.cuda.is_available())
print('GPU:', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'NONE')
if not torch.cuda.is_available():
    raise SystemExit('CUDA not available')
PY

echo "=== 20-STEP PREFLIGHT ==="
rfd3 design \
  out_dir="$PRE" \
  inputs="$SPEC" \
  n_batches=1 \
  diffusion_batch_size=1 \
  inference_sampler.num_timesteps=20 \
  inference_sampler.step_scale=3 \
  inference_sampler.gamma_0=0.2 \
  prevalidate_inputs=true \
  skip_existing=true \
  dump_trajectories=false

PRECOUNT=$(find "$PRE" -type f -name '*.cif.gz' | wc -l)
echo "preflight designs=$PRECOUNT"
if [ "$PRECOUNT" -lt 3 ]; then
  echo "Preflight did not produce all 3 length arms."
  exit 30
fi

echo "============================================================"
echo "PRODUCTION: 3 ARMS x 16 BATCHES x 4 = 192 NEW BACKBONES"
echo "============================================================"
rfd3 design \
  out_dir="$PROD" \
  inputs="$SPEC" \
  n_batches=16 \
  diffusion_batch_size=4 \
  inference_sampler.num_timesteps=200 \
  inference_sampler.step_scale=3 \
  inference_sampler.gamma_0=0.2 \
  low_memory_mode=true \
  prevalidate_inputs=true \
  skip_existing=true \
  dump_trajectories=false

COUNT=$(find "$PROD" -type f -name '*.cif.gz' | wc -l)
echo "generated=$COUNT"

deactivate || true
source "$GEOMENV/bin/activate"
python pipeline/lib/screen_vertical_rfd3.py \
  --root "$PROD" \
  --ring "$SLYB/inputs/target_ring_C11.pdb" \
  --geometry "$GEO" \
  --out "$SCREEN"

echo "=== TOP 25 ==="
column -t -s $'\t' "$SCREEN/screen.tsv" | head -26 || true
NREVIEW=$(awk -F '\t' 'NR>1 && $2=="NEEDS_VISUAL_REVIEW" {n++} END {print n+0}' "$SCREEN/screen.tsv")
echo "generated=$COUNT"
echo "NEEDS_VISUAL_REVIEW=$NREVIEW"

LATEST="$REPO/reports/slyb/rfd3_vertical_latest"
rm -rf "$LATEST"
mkdir -p "$LATEST/top_structures"
cp "$SPEC" "$LATEST/vertical_spec.json"
cp "$GEO" "$LATEST/vertical_geometry.json"
cp "$SCREEN/screen.tsv" "$LATEST/screen.tsv"
find "$SCREEN/top_structures" -type f -name '*.pdb' | head -60 | while IFS= read -r f; do
  cp "$f" "$LATEST/top_structures/"
done

cat > "$LATEST/README.txt" <<EOF
SlyB vertical-groove RFD3 campaign
Job: $SLURM_JOB_ID
Generated new backbones: $COUNT
Candidates needing visual review: $NREVIEW

Goal:
- binder in groove between adjacent native SlyB subunits
- COM guided to recovered design2 groove position
- four interface hotspots, two on each subunit, vertically separated
- binder principal CA axis <=20 degrees from native C11 symmetry axis
- target CA RMSD <=2.0 A
- placement residual <=1.0 A
- >=3 distinct binder residues contacting each partner at 5 A on every copy
- binder-ring backbone minimum >=1.8 A
- all 55 binder-binder pairs minimum >=1.8 A

Raw campaign: $RUN
Screen: $SCREEN/screen.tsv
PyMOL-ready structures: $SCREEN/top_structures
EOF

touch "$RUN/DONE"
echo "============================================================"
echo "FINISHED"
echo "END=$(date)"
echo "RUN=$RUN"
echo "SCREEN=$SCREEN/screen.tsv"
echo "TOP=$SCREEN/top_structures"
echo "============================================================"
