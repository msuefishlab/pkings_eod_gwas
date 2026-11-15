#!/bin/bash
# w_joint.sh — array worker: two-population (joint) vcf2smc for smc++ split.
# Manifest line (TAB): PAIR POP1 POP2 CHR DI_TAG DI_SAMPLE
# 2 DIs x 25 chr per pair; masks = union of the two pops' masks_v2 regions; APA_2009 is
# labelled "APA" so it matches the marginal model pid. Submitted by orchestrate.sh,
# which exports SMCPP_WORK.
set -euo pipefail
MANIFEST=$1; V2=$2; REPO=$3   # passed by orchestrate.sh (array jobs run from a spool copy)
W=${SMCPP_WORK:?run through orchestrate.sh}
INP=$REPO/input_data/12_Historical_Ne
SIF=$REPO/images/smcpp.sif
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
[ -n "$LINE" ] || exit 0
IFS=$'\t' read -r PAIR POP1 POP2 CHR DI_TAG DI_SAMPLE <<< "$LINE"
OUT=$V2/smc_input_split/$PAIR/${CHR}_${DI_TAG}.smc.gz
mkdir -p "$(dirname "$OUT")" "$V2/logs/joint"
if [ -s "$OUT" ]; then echo "[skip] $OUT"; exit 0; fi
lab() { [ "$1" = "APA_2009" ] && echo APA || echo "$1"; }
S1=$(tr '\n' ',' < "$INP/samples/${POP1}_sub20.txt" | sed 's/,$//')
S2=$(tr '\n' ',' < "$INP/samples/${POP2}_sub20.txt" | sed 's/,$//')
VCF=$W/data/vcf_by_chr/${CHR}.vcf.gz
MASK=$V2/masks_split/${POP1}_${POP2}_union.mask.bed.gz
for f in "$VCF" "$MASK"; do [ -s "$f" ] || { echo "[fail] missing $f"; exit 1; }; done
TMP=$(dirname "$OUT")/.tmp.$$.$(basename "$OUT")
t0=$(date +%s)
singularity exec "$SIF" smc++ vcf2smc --mask "$MASK" -d "$DI_SAMPLE" "$DI_SAMPLE" \
    "$VCF" "$TMP" "$CHR" "$(lab $POP1):${S1}" "$(lab $POP2):${S2}" \
    > "$V2/logs/joint/${PAIR}_${CHR}_${DI_TAG}.log" 2>&1
if [ -s "$TMP" ] && gzip -t "$TMP" 2>/dev/null && [ "$(zcat "$TMP" | head -2 | wc -l)" -ge 2 ]; then
  mv -f "$TMP" "$OUT"; echo "[ok] $PAIR $CHR $DI_TAG $(( $(date +%s)-t0 ))s"
else
  rm -f "$TMP"; echo "[fail] $PAIR $CHR $DI_TAG"; exit 1
fi
