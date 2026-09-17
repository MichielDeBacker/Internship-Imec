#!/bin/bash -l
#SBATCH --job-name=rfd3_setup312
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/rfd3_setup_%j.out
#SBATCH --error=/scratch/brussel/vo/000/bvo00014/vsc39230/slyb/logs/rfd3_setup_%j.out

set -Eeuo pipefail

BASE=/scratch/brussel/vo/000/bvo00014/vsc39230
RFDENV="$BASE/venvs/rfd3-foundry"
CKPT="$BASE/models/foundry_checkpoints"

echo "============================================================"
echo "RFD3 SETUP"
echo "JOB=$SLURM_JOB_ID"
echo "HOST=$(hostname)"
echo "START=$(date)"
echo "============================================================"

module purge
module load Python/3.12.3-GCCcore-13.3.0

echo
echo "=== SYSTEM PYTHON ==="
which python
python --version

PYVER=$(python -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')

if [ "$PYVER" != "3.12" ]; then
    echo "ERROR: expected Python 3.12, got $PYVER"
    exit 20
fi

mkdir -p "$BASE/venvs" "$CKPT"

echo
echo "=== HANDLE OLD BROKEN VENV ==="

if [ -d "$RFDENV" ]; then
    OLDVER=$(
        "$RFDENV/bin/python" --version 2>&1 \
        || true
    )

    echo "existing env: $OLDVER"

    if [[ "$OLDVER" != Python\ 3.12* ]]; then
        BACKUP="${RFDENV}.old_$(date +%Y%m%d_%H%M%S)"
        echo "moving old env -> $BACKUP"
        mv "$RFDENV" "$BACKUP"
    fi
fi

echo
echo "=== CREATE / UPDATE PYTHON 3.12 VENV ==="

if [ ! -x "$RFDENV/bin/python" ]; then
    python -m venv "$RFDENV"
fi

"$RFDENV/bin/python" --version

"$RFDENV/bin/python" -m pip install \
    --upgrade \
    pip \
    setuptools \
    wheel

echo
echo "=== INSTALL RFD3 ==="

"$RFDENV/bin/python" -m pip install \
    "rc-foundry[rfd3]"

echo
echo "=== VERIFY EXECUTABLES ==="

test -x "$RFDENV/bin/rfd3"
test -x "$RFDENV/bin/foundry"

"$RFDENV/bin/rfd3" --help >/dev/null
"$RFDENV/bin/foundry" --help >/dev/null

echo "rfd3=$RFDENV/bin/rfd3"
echo "foundry=$RFDENV/bin/foundry"

echo
echo "=== INSTALL CHECKPOINT ==="

export FOUNDRY_CHECKPOINT_DIRS="$CKPT"

"$RFDENV/bin/foundry" install rfd3 \
    --checkpoint-dir "$CKPT"

echo
echo "=== PACKAGE VERSION ==="

"$RFDENV/bin/python" - <<'PY'
import sys
import importlib.metadata

print("python:", sys.version)

for pkg in ("rc-foundry", "torch"):
    try:
        print(pkg, importlib.metadata.version(pkg))
    except Exception as e:
        print(pkg, e)

import torch
print("torch:", torch.__version__)
print("CUDA build:", torch.version.cuda)
PY

"$RFDENV/bin/python" -m pip freeze \
    > "$RFDENV/requirements-lock.txt"

touch "$RFDENV/RFD3_READY"

echo
echo "============================================================"
echo "RFD3 SETUP COMPLETE"
echo "READY=$RFDENV/RFD3_READY"
echo "END=$(date)"
echo "============================================================"
