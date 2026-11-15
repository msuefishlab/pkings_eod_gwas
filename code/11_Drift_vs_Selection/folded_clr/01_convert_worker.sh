#!/bin/bash
# $1 = "POP CHROM"  -- folded SweepFinder2 input for one lineage x chromosome.
# Paths come from common.sh.
set -euo pipefail
read -r pop chrom <<< "$1"
vcf="${SD}/by_chrom/${pop}.polarized.${chrom}.vcf.gz"
out="${SF2DIR}/sfs_input_folded/${pop}.${chrom}.sf2.txt"
[[ -f "$vcf" ]] || { echo "MISSING $vcf"; exit 1; }
if [[ -s "$out" ]]; then echo "skip $pop $chrom"; exit 0; fi
singularity exec "${gwas_tools_image}" \
  python3 "${root}/code/11_Drift_vs_Selection/sweep_detection/sweepfinder2/vcf_to_sweepfinder2.py" \
  "$vcf" "${out}.tmp" -c "$chrom" --polymorphic-only --folded > "$LOGDIR/convert.${pop}.${chrom}.log" 2>&1
mv "${out}.tmp" "$out"
echo "done $pop $chrom $(($(wc -l < "$out") - 1)) sites"
