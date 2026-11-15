#!/bin/bash
# $1 = POP  -- concatenate 25 chromosomes and compute the folded background spectrum.
# The unfolded spectrum is NOT reused (METHODS A4).
set -euo pipefail
pop="$1"
indir="${SF2DIR}/sfs_input_folded"
outdir="${SF2DIR}/background_folded"
combined="$outdir/${pop}.allchroms.sf2.txt"
spectrum="$outdir/${pop}.spectrum"
head -1 "$indir/${pop}.chr1.sf2.txt" > "$combined"
for i in $(seq 1 25); do
  f="$indir/${pop}.chr${i}.sf2.txt"
  [[ -f "$f" ]] || { echo "ERROR missing $f"; exit 1; }
  tail -n +2 "$f" >> "$combined"
done
echo "$pop combined sites: $(($(wc -l < "$combined") - 1))"
singularity exec "${sweepfinder2_image}" \
  SweepFinder2 -f "$combined" "${spectrum}.tmp" > "$LOGDIR/background.${pop}.log" 2>&1
mv "${spectrum}.tmp" "$spectrum"
echo "$pop spectrum lines: $(wc -l < "$spectrum")"
