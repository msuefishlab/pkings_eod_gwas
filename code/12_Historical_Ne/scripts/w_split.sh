#!/bin/bash
# w_split.sh — array worker for smc++ split.
# Manifest line (TAB): PAIR POP1 POP2 MU_LABEL MODE   (MODE = unfolded | folded)
# unfolded = no polarization flag;
# folded   = adds --polarization-error 0.5. smc++ 1.15.4 split ignores this flag
#            (identical starting log-likelihoods), so the two modes are replicate fits.
set -euo pipefail
MANIFEST=$1; V2=$2; REPO=$3   # passed by orchestrate.sh (array jobs run from a spool copy)
SIF=$REPO/images/smcpp.sif
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
[ -n "$LINE" ] || exit 0
IFS=$'\t' read -r PAIR POP1 POP2 MU MODE <<< "$LINE"
SUB=split; EXTRA=""
[ "$MODE" = folded ] && { SUB=split_folded; EXTRA="--polarization-error 0.5"; }
OUT=$V2/$SUB/${PAIR}_${MU}
if [ -s "$OUT/model.final.json" ]; then echo "[skip] $OUT"; exit 0; fi
M1=$V2/estimates/${POP1}_${MU}/model.final.json
M2=$V2/estimates/${POP2}_${MU}/model.final.json
for f in "$M1" "$M2"; do [ -s "$f" ] || { echo "[fail] missing $f"; exit 1; }; done
N=$(ls "$V2/smc_input_split/$PAIR"/*.smc.gz 2>/dev/null | wc -l)
[ "$N" -eq 50 ] || { echo "[fail] $PAIR joint inputs $N/50"; exit 1; }
rm -rf "$OUT"; mkdir -p "$OUT"
singularity exec "$SIF" smc++ split --cores ${SLURM_CPUS_PER_TASK:-8} $EXTRA -o "$OUT" "$M1" "$M2" \
    "$V2/smc_input_split/$PAIR"/*.smc.gz > "$OUT/split.log" 2>&1
[ -s "$OUT/model.final.json" ] || { echo "[fail] no split model $OUT"; exit 1; }
echo "[ok] $PAIR $MU $MODE"
