#!/bin/bash
# w_estimate.sh — array worker for point estimates and bootstraps.
# Manifest line (TAB): KIND K POP MU MU_LABEL     (KIND = point | boot)
set -euo pipefail
MANIFEST=$1; V2=$2; REPO=$3   # passed by orchestrate.sh (array jobs run from a spool copy)
SIF=$REPO/images/smcpp.sif
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
[ -n "$LINE" ] || exit 0
IFS=$'\t' read -r KIND K POP MU LABEL <<< "$LINE"
if [ "$KIND" = point ]; then
  IN=$V2/smc_input/$POP; OUT=$V2/estimates/${POP}_mu${LABEL}
else
  IN=$V2/bootstraps/boot_${K}/${POP}_mu${LABEL}; OUT=$V2/bootstraps/estimates/boot_${K}/${POP}_mu${LABEL}
fi
if [ -s "$OUT/model.final.json" ]; then echo "[skip] $OUT"; exit 0; fi
N=$(ls "$IN"/chr*.smc.gz 2>/dev/null | wc -l)
[ "$N" -eq 125 ] || { echo "[fail] $IN has $N/125 inputs"; exit 1; }
rm -rf "$OUT"; mkdir -p "$OUT"
echo "$(date): $KIND K=$K POP=$POP MU=$MU host=$(hostname) job=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}" > "$OUT/estimate.log"
singularity run "$SIF" estimate --cores 4 --polarization-error 0.5 -o "$OUT" "$MU" "$IN"/chr*.smc.gz >> "$OUT/estimate.log" 2>&1
[ -s "$OUT/model.final.json" ] || { echo "[fail] no model.final.json in $OUT"; exit 1; }
echo "$(date): DONE" >> "$OUT/estimate.log"; echo "[ok] $KIND $K $POP $LABEL"
