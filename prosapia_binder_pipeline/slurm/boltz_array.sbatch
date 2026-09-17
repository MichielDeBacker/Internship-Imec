#!/bin/bash
#SBATCH --job-name=slyb_boltz
#SBATCH --partition=hopper_gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=logs/boltz_%A_%a.out
#SBATCH --error=logs/boltz_%A_%a.err

set -euo pipefail

PIPELINE_ROOT="${SLURM_SUBMIT_DIR:-/scratch/brussel/vo/000/bvo00014/vsc39230/repos/Internship-Imec/prosapia_binder_pipeline}"
CONFIG="${CONFIG:-$PIPELINE_ROOT/config/slyb_hydra.yaml}"
PIPELINE_PYTHON="${PIPELINE_PYTHON:-/data/leuven/392/vsc39230/slyb_envs/prosapia-pipeline/bin/python}"

VSC_DATA=/data/leuven/392/vsc39230
CACHE="$VSC_DATA/slyb_prosapia_cache/boltz_runtime/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$CACHE/tmp" "$CACHE/xdg" "$CACHE/triton" "$CACHE/cuda"
export TMPDIR="$CACHE/tmp"
export XDG_CACHE_HOME="$CACHE/xdg"
export TRITON_CACHE_DIR="$CACHE/triton"
export CUDA_CACHE_PATH="$CACHE/cuda"

cd "$PIPELINE_ROOT"
"$PIPELINE_PYTHON" scripts/05_boltz_one.py \
    --config "$CONFIG" \
    --index "$SLURM_ARRAY_TASK_ID"
