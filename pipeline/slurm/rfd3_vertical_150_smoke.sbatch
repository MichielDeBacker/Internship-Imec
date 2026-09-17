#!/bin/bash -l
#SBATCH --job-name=slyb_v150_smoke
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=12
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/v150_smoke_%j.out
#SBATCH --error=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/v150_smoke_%j.out

set -Eeuo pipefail

BASE=/scratch/brussel/vo/000/bvo00014/vsc39230
SLYB="$BASE/slyb"
REPO="$BASE/repos/Internship-Imec"
RFDENV="$BASE/venvs/rfd3-foundry"
CKPT="$BASE/models/foundry_checkpoints"
WORK="$REPO/work/vertical_150_preflight"

RUN="$SLYB/out/rfd3_vertical150_smoke_${SLURM_JOB_ID}"

module purge
module load Python/3.12.3-GCCcore-13.3.0

if [ ! -f "$RFDENV/RFD3_READY" ]; then
    echo "RFD3 environment not ready"
    exit 10
fi

if [ ! -f "$WORK/CPU_PREFLIGHT_OK" ]; then
    echo "CPU geometry preflight missing"
    exit 11
fi

source "$RFDENV/bin/activate"

export FOUNDRY_CHECKPOINT_DIRS="$CKPT"

mkdir -p "$RUN"
cd "$REPO"

echo "============================================================"
echo "RFD3 SLYB SMOKE"
echo "JOB=$SLURM_JOB_ID"
echo "HOST=$(hostname)"
echo "START=$(date)"
echo "============================================================"

nvidia-smi

python - <<'PY'
import torch

print("torch:", torch.__version__)
print("cuda available:", torch.cuda.is_available())

if not torch.cuda.is_available():
    raise SystemExit("CUDA unavailable")

print("GPU:", torch.cuda.get_device_name(0))
print(
    "VRAM GB:",
    torch.cuda.get_device_properties(0).total_memory / 1024**3
)
PY

echo
echo "=== PREVALIDATE + GENERATE 3 BACKBONES ==="

rfd3 design \
    inputs="$WORK/vertical_spec.json" \
    out_dir="$RUN" \
    n_batches=1 \
    diffusion_batch_size=1 \
    inference_sampler.num_timesteps=20 \
    inference_sampler.step_scale=3 \
    inference_sampler.gamma_0=0.2 \
    prevalidate_inputs=True \
    low_memory_mode=False \
    skip_existing=False \
    dump_trajectories=False

COUNT=$(
    find "$RUN" \
        -type f \
        -name '*_model_0.cif.gz' \
        | wc -l
)

echo
echo "smoke generated=$COUNT"
echo "smoke expected=3"

if [ "$COUNT" -ne 3 ]; then
    echo "ERROR: smoke did not produce exactly 3 backbones"
    exit 30
fi

touch "$RUN/SMOKE_OK"

echo
echo "SMOKE PASS"
echo "END=$(date)"
