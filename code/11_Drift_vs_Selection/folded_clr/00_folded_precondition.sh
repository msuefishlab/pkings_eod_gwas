#!/bin/bash
## 00_folded_precondition.sh -- METHODS A3: verify the folded conversion is genuinely
## polarization-free before the genome-wide scan is run, i.e. that it does not treat REF as
## ancestral.
##
## Checks P1-P5 on chr1 for lineage BP2. Exits non-zero if any required check fails, which
## stops run_all.sh before a single scan is launched.
set -euo pipefail
cd "$(dirname "$0")"; source common.sh
POP=${1:-BP2}; CHROM=${2:-chr1}

vcf="${SD}/by_chrom/${POP}.polarized.${CHROM}.vcf.gz"
[ -s "$vcf" ] || { echo "MISSING ${vcf} -- run code/11_Drift_vs_Selection/sweep_detection/ehh/01b_split_by_chrom_array.sb first" >&2; exit 1; }

## The check compares the folded input against the unfolded one site for site (P2, P3), so both
## must exist for this chromosome. Both are regenerated here if absent, with the repository's
## own converter and the same flags as sweepfinder2/10_sf2_convert.sb, so the precondition is executable
## from a clean clone.
conv="${root}/code/11_Drift_vs_Selection/sweep_detection/sweepfinder2/vcf_to_sweepfinder2.py"
for mode in unfolded folded; do
  out="${SF2DIR}/sfs_input_${mode}/${POP}.${CHROM}.sf2.txt"
  mkdir -p "$(dirname "$out")"
  if [ -s "$out" ]; then echo "have ${out}"; continue; fi
  extra=""; [ "$mode" = folded ] && extra="--folded"
  echo "building ${out}"
  sx python3 "$conv" "$vcf" "${out}.tmp" -c "${CHROM}" --polymorphic-only ${extra} \
      > "${LOGDIR}/precondition_convert.${POP}.${CHROM}.${mode}.log" 2>&1
  mv "${out}.tmp" "$out"
done

sx python3 check_folded_conversion.py "${POP}" "${CHROM}" "${FCLR_OUT}/audit/folded_conversion_check.tsv"
cp -f "${FCLR_OUT}/audit/folded_conversion_check.tsv" "${FCLR_OUT}/folded_conversion_check.tsv"
