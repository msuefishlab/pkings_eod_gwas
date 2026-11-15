#!/bin/bash
# Usage: bash w_vcf2smc.sh VCF_GZ OUTFILE CHR DI_SAMPLE SMCSIF SAMPLE_LIST
VCF=$1; OUTF=$2; CHR=$3; DI=$4; SMCSIF=$5; SAMPLE_LIST=$6
if [ -f "$OUTF" ]; then
    echo "skip $CHR $DI"
    exit 0
fi
/usr/bin/singularity exec "$SMCSIF" smc++ vcf2smc \
    --cores 1 \
    -d "$DI" "$DI" \
    "$VCF" "$OUTF" "$CHR" \
    "POP:$SAMPLE_LIST" 2>&1 | tail -3
RC=$?
if [ $RC -ne 0 ]; then
    echo "ERROR ($RC) in vcf2smc for $CHR $DI" >&2
    exit 0
fi
echo "done $CHR $DI"
