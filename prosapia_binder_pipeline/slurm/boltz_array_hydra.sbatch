#!/bin/bash
#SBATCH --job-name=slyb_boltz_hydra_prod
#SBATCH --partition=hopper_gpu,ampere_gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/brussel/vo/000/bvo00014/vsc39230/repos/Internship-Imec/prosapia_binder_pipeline/logs/hydra_prod_%A_%a.out
#SBATCH --error=/scratch/brussel/vo/000/bvo00014/vsc39230/repos/Internship-Imec/prosapia_binder_pipeline/logs/hydra_prod_%A_%a.err
#SBATCH --chdir=/scratch/brussel/vo/000/bvo00014/vsc39230/repos/Internship-Imec/prosapia_binder_pipeline

module purge
module use /apps/brussel/RL9/zen2-ib/modules/system/all
module use /apps/brussel/RL9/zen2-ib/modules/2024a/all
module load Boltz-2/2.2.1-foss-2024a-CUDA-12.6.0

export TMPDIR=/tmp
export TMP=/tmp
export TEMP=/tmp

CACHE="/data/leuven/392/vsc39230/slyb_prosapia_cache/boltz_runtime/hydra_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export XDG_CACHE_HOME="$CACHE/xdg"
export NUMBA_CACHE_DIR="$CACHE/numba"
export MPLCONFIGDIR="$CACHE/matplotlib"
export TRITON_CACHE_DIR="$CACHE/triton"
export CUDA_CACHE_PATH="$CACHE/cuda"

mkdir -p \
    "$XDG_CACHE_HOME" \
    "$NUMBA_CACHE_DIR" \
    "$MPLCONFIGDIR" \
    "$TRITON_CACHE_DIR" \
    "$CUDA_CACHE_PATH"

/data/leuven/392/vsc39230/slyb_envs/prosapia-pipeline/bin/python \
    scripts/05_boltz_one.py \
    --config config/slyb_hydra.yaml \
    --index "$SLURM_ARRAY_TASK_ID"
