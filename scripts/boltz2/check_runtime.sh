#!/usr/bin/env bash
set -euo pipefail

eval "$(/software/AMD/Conda/miniforge3/bin/conda shell.bash hook)"
conda activate boltz2

export BOLTZ_CACHE=/software/General-Models/Boltz2/

echo "Python: $(command -v python)"
echo "Boltz:  $(command -v boltz)"
echo "Cache:  $BOLTZ_CACHE"

python - <<'PY'
import importlib.metadata
import torch

print("Boltz:", importlib.metadata.version("boltz"))
print("PyTorch:", torch.__version__)
print("CUDA:", torch.cuda.is_available())
PY

nvidia-smi \
  --query-gpu=index,memory.used,memory.total,utilization.gpu \
  --format=csv,noheader
