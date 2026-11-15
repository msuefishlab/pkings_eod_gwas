#!/bin/bash
# w_vcf2smc.sh — array worker: single-population vcf2smc on masks_v2.
# Manifest line (TAB): POP CHR DI DIIDX
# Output tree = $V2/smc_input. Submitted by orchestrate.sh, which exports SMCPP_WORK.
set -euo pipefail
MANIFEST=$1; V2=$2; REPO=$3   # passed by orchestrate.sh (array jobs run from a spool copy)
W=${SMCPP_WORK:?run through orchestrate.sh}
INP=$REPO/input_data/12_Historical_Ne
SIF=$REPO/images/smcpp.sif
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
[ -n "$LINE" ] || { echo "empty line ${SLURM_ARRAY_TASK_ID}"; exit 0; }
IFS=$'\t' read -r POP CHR DI DIIDX <<< "$LINE"
OUT=$V2/smc_input/$POP/${CHR}_di${DIIDX}.smc.gz
mkdir -p "$(dirname "$OUT")" "$V2/logs/vcf2smc"
if [ -s "$OUT" ]; then echo "[skip] $OUT"; exit 0; fi
SMP_POP=$POP; [ "$POP" = "APA_2009" ] && SMP_POP=APA
SAMPLES=$(tr '\n' ',' < "$INP/samples/${POP}_sub20.txt" | sed 's/,$//')
VCF=$W/data/vcf_by_chr/${CHR}.vcf.gz
MASK=$W/audit/masks_v2/${POP}.mask.bed.gz
for f in "$VCF" "$MASK"; do [ -s "$f" ] || { echo "[fail] missing $f"; exit 1; }; done
TMP=$(dirname "$OUT")/.tmp.$$.$(basename "$OUT")      # keep .smc.gz suffix (gzip by extension)
t0=$(date +%s)
singularity exec "$SIF" smc++ vcf2smc --mask "$MASK" -d "$DI" "$DI" \
    "$VCF" "$TMP" "$CHR" "${SMP_POP}:${SAMPLES}" \
    > "$V2/logs/vcf2smc/${POP}_${CHR}_di${DIIDX}.log" 2>&1
if [ -s "$TMP" ] && gzip -t "$TMP" 2>/dev/null && [ "$(zcat "$TMP" | head -2 | wc -l)" -ge 2 ]; then
  mv -f "$TMP" "$OUT"; echo "[ok] $POP $CHR di$DIIDX $(( $(date +%s)-t0 ))s $(stat -c %s "$OUT")B"
else
  rm -f "$TMP"; echo "[fail] $POP $CHR di$DIIDX empty/corrupt output"; exit 1
fi
