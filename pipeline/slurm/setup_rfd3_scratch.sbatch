#!/bin/bash -l
#SBATCH --job-name=rfd3_setup312
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/rfd3_setup_%j.out
#SBATCH --error=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/rfd3_setup_%j.out

set -Eeuo pipefail

BASE=/scratch/brussel/vo/000/bvo00014/vsc39230

RFDENV="$BASE/venvs/rfd3-foundry"
CKPT="$BASE/models/foundry_checkpoints"

export TMPDIR="$BASE/tmp/rfd3_tmp"
export XDG_CACHE_HOME="$BASE/tmp/rfd3_xdg"
export PIP_CACHE_DIR="$BASE/tmp/rfd3_pip"
export PIP_NO_CACHE_DIR=1

mkdir -p \
    "$TMPDIR" \
    "$XDG_CACHE_HOME" \
    "$PIP_CACHE_DIR" \
    "$CKPT" \
    "$BASE/venvs"

module purge
module load Python/3.12.3-GCCcore-13.3.0

echo "system python:"
python --version
which python

if [ -x "$RFDENV/bin/python" ]; then

    OLD=$(
        "$RFDENV/bin/python" \
        --version \
        2>&1 \
        || true
    )

    echo "existing RFDENV: $OLD"

    if [[ "$OLD" != Python\ 3.12* ]]
    then
        mv \
            "$RFDENV" \
            "${RFDENV}.old_$(date +%Y%m%d_%H%M%S)"
    fi
fi

if [ ! -x "$RFDENV/bin/python" ]; then
    python -m venv "$RFDENV"
fi

"$RFDENV/bin/python" \
    -m pip install \
    --upgrade \
    pip setuptools wheel

"$RFDENV/bin/python" \
    -m pip install \
    "rc-foundry[rfd3]"

test -x "$RFDENV/bin/rfd3"
test -x "$RFDENV/bin/foundry"

export FOUNDRY_CHECKPOINT_DIRS="$CKPT"

"$RFDENV/bin/foundry" \
    install rfd3 \
    --checkpoint-dir "$CKPT"

"$RFDENV/bin/rfd3" \
    --help >/dev/null

touch "$RFDENV/RFD3_READY"

echo
echo "RFD3 READY"
echo "python:"
"$RFDENV/bin/python" --version
echo "rfd3:"
"$RFDENV/bin/rfd3" --help \
    | head -20
