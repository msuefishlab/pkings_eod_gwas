#!/bin/bash
# $1 = "POP CHROM"  -- SweepFinder2 folded CLR scan for one lineage x chromosome, -lg 10000
# (the 10-kb grid of the unfolded scan, sweepfinder2/12_sf2_scan.sb).
set -euo pipefail
read -r pop chrom <<< "$1"
input="${SF2DIR}/sfs_input_folded/${pop}.${chrom}.sf2.txt"
spectrum="${SF2DIR}/background_folded/${pop}.spectrum"
out="${SF2DIR}/clr_folded/${pop}.${chrom}.clr.txt"
[[ -f "$input"    ]] || { echo "MISSING $input"; exit 1; }
[[ -f "$spectrum" ]] || { echo "MISSING $spectrum"; exit 1; }
if [[ -s "$out" ]]; then echo "skip $pop $chrom"; exit 0; fi
singularity exec "${sweepfinder2_image}" \
  SweepFinder2 -lg 10000 "$input" "$spectrum" "${out}.tmp" > "$LOGDIR/scan.${pop}.${chrom}.log" 2>&1
mv "${out}.tmp" "$out"
echo "done $pop $chrom $(wc -l < "$out") grid points"
