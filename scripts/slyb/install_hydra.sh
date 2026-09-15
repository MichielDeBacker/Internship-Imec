#!/bin/bash
# Install RFdiffusion3 + LigandMPNN on hydra (VUB Tier-2, Rocky Linux 9, Lmod).
#
# WHY MODULES AND NOT A CONTAINER -- decided by probe, not preference:
#   1. rc-foundry requires Python >=3.12,<3.13. The module
#      PyTorch/2.6.0-foss-2024a-CUDA-12.6.0 provides exactly Python 3.12.3
#      with torch 2.6.0 / CUDA 12.6 -- inside that window. The newer
#      foss-2025a stack ships Python 3.13.1, which rc-foundry EXCLUDES.
#   2. Apptainer 1.5.3 exists but the user has no /etc/subuid entry, so
#      `apptainer build --fakeroot` is unavailable. A container would have to
#      be built off-cluster and a ~8 GB image shipped in, with no upside.
#   3. A venv with --system-site-packages reuses the module's CUDA-linked
#      torch, so only pure-Python design code is installed on top (~400 MB
#      instead of ~5 GB) and the vendor-tuned torch build is preserved.
#
# Run on the LOGIN node (compute nodes have no outbound network).
set -euo pipefail

ROOT="${SLYB_ROOT:-/scratch/brussel/vo/000/bvo00014/vsc39230/slyb}"
VENV="$ROOT/venv-rfd3"
CKPT="$ROOT/rfd3ckpt"
TORCH_MODULE="PyTorch/2.6.0-foss-2024a-CUDA-12.6.0"

echo "=== target layout ==="
echo "  root        $ROOT"
echo "  venv        $VENV"
echo "  checkpoints $CKPT"

mkdir -p "$ROOT" "$CKPT" "$ROOT/inputs" "$ROOT/out" "$ROOT/logs"

source /etc/profile
module purge
module load "$TORCH_MODULE"

PYV=$(python -c 'import sys;print("%d.%d"%sys.version_info[:2])')
echo "=== module python $PYV ==="
python - <<'PYCHK'
import sys
if not ((3,12) <= sys.version_info[:2] < (3,13)):
    sys.exit("FATAL: rc-foundry needs Python >=3.12,<3.13; module gave %d.%d"
             % sys.version_info[:2])
print("python version inside rc-foundry window: OK")
PYCHK

# tmpfs check -- checkpoints on RAM-backed storage get the run OOM-killed
TMPTYPE=$(df -T /tmp | tail -1 | awk '{print $2}')
echo "=== /tmp is $TMPTYPE ==="
[ "$TMPTYPE" = "tmpfs" ] && { echo "FATAL: /tmp is tmpfs"; exit 1; }

if [ ! -d "$VENV" ]; then
    echo "=== creating venv on top of module torch ==="
    python -m venv --system-site-packages "$VENV"
fi
# shellcheck disable=SC1091
source "$VENV/bin/activate"
python -m pip install --upgrade pip wheel

echo "=== verifying module torch is visible inside venv ==="
python - <<'PYCHK'
import torch, sys
print("torch", torch.__version__, "cuda-built", torch.version.cuda)
if torch.version.cuda is None:
    sys.exit("FATAL: venv resolved a CPU-only torch; --system-site-packages lost")
PYCHK

echo "=== installing rc-foundry[rfd3] ==="
# VERIFIED by dry-run on hydra: resolving rc-foundry[rfd3] against the module
# torch does NOT list torch among the packages to install -- pip accepts the
# module's torch 2.6.0 as satisfying "torch<3,>=2.2.0". So no pin is needed and
# none is applied; forcing one would risk shadowing the CUDA-linked build.
# 74 pure-Python packages are added (atomworks, biotite, rdkit, lightning, ...).
python -m pip install --no-cache-dir "rc-foundry[rfd3]"

echo "=== confirming pip did not replace torch ==="
python - <<'PYCHK'
import torch, sys
print("torch after install:", torch.__version__, "cuda-built", torch.version.cuda)
if torch.version.cuda is None:
    sys.exit("FATAL: rc-foundry install pulled a CPU-only torch wheel")
PYCHK

echo "=== downloading RFD3 checkpoints (~2.6 GB) ==="
export FOUNDRY_CHECKPOINT_DIRS="$CKPT"
if [ -f "$CKPT/rfd3_latest.ckpt" ]; then
    echo "  checkpoint already present, skipping"
else
    foundry install rfd3 --checkpoint-dir "$CKPT"
fi

echo "=== installing LigandMPNN ==="
LMP="$ROOT/LigandMPNN"
if [ ! -d "$LMP" ]; then
    git clone https://github.com/dauparas/LigandMPNN.git "$LMP"
fi
cd "$LMP"
python -m pip install --no-cache-dir -r requirements.txt || \
    echo "WARN: LigandMPNN requirements partially failed -- check manually"
if [ ! -f "$LMP/model_params/ligandmpnn_v_32_010_25.pt" ]; then
    bash get_model_params.sh "$LMP/model_params"
fi

cat > "$ROOT/env.sh" <<EOF
# source this before any RFD3/LigandMPNN work on hydra
source /etc/profile
module purge
module load $TORCH_MODULE
source $VENV/bin/activate
export FOUNDRY_CHECKPOINT_DIRS=$CKPT
export SLYB_ROOT=$ROOT
export LIGANDMPNN_DIR=$LMP
EOF

echo
echo "=== install complete ==="
echo "  source $ROOT/env.sh   then run verify_hydra.py"
