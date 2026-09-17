#!/bin/bash -l
#SBATCH --job-name=slyb_v150
#SBATCH --array=0-5%6
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=12
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/v150_%A_%a.out
#SBATCH --error=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/v150_%A_%a.out

set -Eeuo pipefail

BASE=/scratch/brussel/vo/000/bvo00014/vsc39230
SLYB="$BASE/slyb"
REPO="$BASE/repos/Internship-Imec"
RFDENV="$BASE/venvs/rfd3-foundry"
CKPT="$BASE/models/foundry_checkpoints"
WORK="$REPO/work/vertical_150_preflight"

TASK="$SLURM_ARRAY_TASK_ID"

case "$TASK" in
    0|1)
        KEY="vertical_46_54"
        ;;
    2|3)
        KEY="vertical_54_62"
        ;;
    4|5)
        KEY="vertical_62_72"
        ;;
    *)
        echo "invalid array task: $TASK"
        exit 10
        ;;
esac

SPEC="$WORK/shards/shard_${TASK}_${KEY}.json"

RUNROOT="$SLYB/out/rfd3_vertical150_${SLURM_ARRAY_JOB_ID}"
OUT="$RUNROOT/task_${TASK}_${KEY}"

module purge
module load Python/3.12.3-GCCcore-13.3.0

test -f "$RFDENV/RFD3_READY"
test -f "$WORK/CPU_PREFLIGHT_OK"
test -s "$SPEC"

source "$RFDENV/bin/activate"

export FOUNDRY_CHECKPOINT_DIRS="$CKPT"

mkdir -p "$OUT"
cd "$REPO"

echo "============================================================"
echo "SLYB RFD3 PRODUCTION"
echo "ARRAY=$SLURM_ARRAY_JOB_ID"
echo "TASK=$TASK"
echo "KEY=$KEY"
echo "HOST=$(hostname)"
echo "START=$(date)"
echo "SPEC=$SPEC"
echo "OUT=$OUT"
echo "============================================================"

nvidia-smi

python - <<'PY'
import torch

if not torch.cuda.is_available():
    raise SystemExit("CUDA unavailable")

print("GPU:", torch.cuda.get_device_name(0))
print(
    "VRAM GB:",
    torch.cuda.get_device_properties(0).total_memory / 1024**3
)
PY

echo
echo "=== GENERATE 25 BACKBONES ON THIS GPU ==="

# 5 batches x 5 diffusion samples = 25
# 6 array tasks x 25 = 150 total
#
# Length allocation:
# tasks 0,1 -> 46-54 = 50
# tasks 2,3 -> 54-62 = 50
# tasks 4,5 -> 62-72 = 50

rfd3 design \
    inputs="$SPEC" \
    out_dir="$OUT" \
    n_batches=5 \
    diffusion_batch_size=5 \
    inference_sampler.num_timesteps=200 \
    inference_sampler.step_scale=3 \
    inference_sampler.gamma_0=0.2 \
    prevalidate_inputs=True \
    low_memory_mode=False \
    skip_existing=False \
    dump_trajectories=False

COUNT=$(
    find "$OUT" \
        -type f \
        -name '*_model_0.cif.gz' \
        | wc -l
)

echo
echo "generated=$COUNT"
echo "expected=25"

if [ "$COUNT" -ne 25 ]; then
    echo "ERROR: shard output count != 25"
    exit 30
fi

touch "$OUT/SHARD_DONE"

echo
echo "SHARD COMPLETE"
echo "END=$(date)"
