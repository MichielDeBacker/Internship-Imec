#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${CONFIG:-$ROOT/config/slyb_hydra.yaml}"
PYTHON="${PIPELINE_PYTHON:-python3}"

usage() {
    cat <<'EOF'
Usage: ./workbench.sh COMMAND

Commands:
  check           Check external paths/tools without modifying Prosapia
  prepare         Register/canonicalize all 150 RFD3 backbones
  submit-seq      Submit ProteinMPNN array for all backbone manifest rows
  collect-seq     Collect ProteinMPNN outputs into sequences.tsv
  prepare-boltz   Build Boltz YAMLs with A/B target templates, C untemplated
  submit-boltz    Submit Boltz array for all prepared sequence inputs
  collect-boltz   Collect ipTM and pairwise ipTM tables
  status          Show campaign stage counts
EOF
}

cmd="${1:-}"
case "$cmd" in
  check)
    "$PYTHON" "$ROOT/scripts/00_check_environment.py" --config "$CONFIG"
    ;;
  prepare)
    "$PYTHON" "$ROOT/scripts/01_prepare_backbones.py" --config "$CONFIG"
    ;;
  submit-seq)
    CONFIG="$CONFIG" bash "$ROOT/scripts/submit_sequence.sh"
    ;;
  collect-seq)
    "$PYTHON" "$ROOT/scripts/03_collect_sequences.py" --config "$CONFIG"
    ;;
  prepare-boltz)
    "$PYTHON" "$ROOT/scripts/04_prepare_boltz_inputs.py" --config "$CONFIG"
    ;;
  submit-boltz)
    CONFIG="$CONFIG" bash "$ROOT/scripts/submit_boltz.sh"
    ;;
  collect-boltz)
    "$PYTHON" "$ROOT/scripts/06_collect_boltz_scores.py" --config "$CONFIG"
    ;;
  status)
    "$PYTHON" "$ROOT/scripts/07_status.py" --config "$CONFIG"
    ;;
  *)
    usage
    exit 2
    ;;
esac
