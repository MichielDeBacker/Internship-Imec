#!/bin/bash
set -uo pipefail
BASE=/scratch/brussel/vo/000/bvo00014/vsc39230
REPO="$BASE/repos/Internship-Imec"
SLYB="$BASE/slyb"
cd "$REPO" || exit 1
mkdir -p "$SLYB/logs"
source "$BASE/venvs/slyb-geometry/bin/activate" || exit 2
python "$REPO/pipeline/lib/build_vertical_groove_spec.py" --help >/dev/null || exit 3
python "$REPO/pipeline/lib/screen_vertical_rfd3.py" --help >/dev/null || exit 4
bash -n "$REPO/pipeline/slurm/rfd3_vertical_overnight.sbatch" || exit 5
SUBMIT=""
for PART in hopper_gpu ampere_gpu; do
  echo "Trying $PART ..."
  SUBMIT=$(sbatch -A bvo00014 --partition="$PART" "$REPO/pipeline/slurm/rfd3_vertical_overnight.sbatch" 2>&1)
  RC=$?
  echo "$SUBMIT"
  if [ "$RC" -eq 0 ] && printf '%s\n' "$SUBMIT" | grep -q 'Submitted batch job'; then
    break
  fi
done
JOBID=$(printf '%s\n' "$SUBMIT" | awk '/Submitted batch job/ {print $4}')
echo "JOBID=$JOBID"
if [ -n "$JOBID" ]; then
  squeue -j "$JOBID" -o '%.18i %.12P %.25j %.10T %.10M %.30R'
  echo "LOG=$SLYB/logs/vertical_rfd3_${JOBID}.out"
else
  echo 'SUBMISSION FAILED'
  exit 6
fi
