#!/bin/bash
# Usage: bash w_estimate.sh SC VALDIR SMCSIF
SC=$1; VALDIR=$2; SMCSIF=$3
SIMDIR="$VALDIR/sim_$SC"
ESTDIR="$VALDIR/estimates/$SC"
mkdir -p "$ESTDIR"
N_SMC=$(ls "$SIMDIR"/*.smc.gz 2>/dev/null | wc -l || echo 0)
if [ "$N_SMC" -eq 0 ]; then
    echo "ERROR: no smc.gz files for $SC" >&2
    exit 0
fi
echo "  $SC: estimating on $N_SMC smc.gz files at $(date)..."
/usr/bin/singularity run "$SMCSIF" estimate \
    --polarization-error 0.5 \
    --cores 4 \
    -o "$ESTDIR/" \
    1e-8 \
    "$SIMDIR"/*.smc.gz 2>&1 | tail -10
echo "  Done: $SC  $(date)"
if [ -f "$ESTDIR/model.final.json" ]; then
    echo "  model.final.json: OK"
else
    echo "  WARNING: model.final.json not found!"
fi
