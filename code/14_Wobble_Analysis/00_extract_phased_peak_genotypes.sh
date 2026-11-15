#!/bin/bash
# Extract phased genotypes at the six GWAS peak SNPs for all focal samples.
# Outputs the phased GT matrix that 01_call_peak_haplotypes.R reads.
#
# Run from repo root after sourcing config.env:
#   source config.env
#   bash code/14_Wobble_Analysis/00_extract_phased_peak_genotypes.sh

set -euo pipefail

source config.env

ANALYSIS="PKINGS_ALL_WOB_EXCLUDED"
SNPS_IN_PEAKS="${root}/output_data/06_Association/${ANALYSIS}/${ANALYSIS}_SNPS_IN_PEAKS.txt"
OUTDIR="${root}/output_data/14_Wobble_Analysis/PKINGS_WOB_MULTILOCUS"
mkdir -p "${OUTDIR}"

# Focal sample list: all BAM + APA + BEN + BIK + BAVA + BIR individuals
# Derived from fish_data_2024a.txt (POP column, col 9)
SAMPLE_LIST="${OUTDIR}/focal_samples.txt"
awk -F'\t' 'NR>1 && $9 ~ /^(BAM|APA|BEN|BIK|BAVA|BIR|DOV|DOG)$/ {print $8}' \
    "${root}/input_data/01_Terra/data_model/fish_data_2024a.txt" \
    > "${SAMPLE_LIST}"
echo "Focal samples written: $(wc -l < "${SAMPLE_LIST}")"

# Build a regions file (tab-delimited CHROM<TAB>POS) from the SNPS_IN_PEAKS file.
# bcftools -R FILE requires tab-delimited columns, not chr:beg-end region strings.
REGIONS_FILE="${OUTDIR}/peak_snp_regions.txt"
awk 'NR>1 {print $1"\t"$3}' "${SNPS_IN_PEAKS}" > "${REGIONS_FILE}"
echo "Peak SNP regions: $(wc -l < "${REGIONS_FILE}")"

# Extract phased VCF subset
PHASED_VCF_OUT="${OUTDIR}/PKINGS_WOB_peak_snps.phased.vcf.gz"
echo "Extracting phased genotypes from ${phased_VCF} ..."

bcftools view \
    -R "${REGIONS_FILE}" \
    -S "${SAMPLE_LIST}" \
    --output-type z \
    -o "${PHASED_VCF_OUT}" \
    "${phased_VCF}"

bcftools index -t "${PHASED_VCF_OUT}"

echo "Phased VCF written: ${PHASED_VCF_OUT}"

# Dump GT matrix (tab-delimited, phased strings preserved: 0|0, 0|1, 1|0, 1|1)
GT_MATRIX="${OUTDIR}/PKINGS_WOB_peak_snps.gt_phased.txt"
bcftools query \
    -H \
    -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
    "${PHASED_VCF_OUT}" \
    > "${GT_MATRIX}"

echo "Phased GT matrix written: ${GT_MATRIX}"
echo "Done. Load ${GT_MATRIX} in R for haplotype-phase analysis."
