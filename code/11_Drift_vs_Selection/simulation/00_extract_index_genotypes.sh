#!/bin/bash
## 00_extract_index_genotypes.sh -- positive control: re-extract the Peak 1/4/5 index-SNP genotypes
## for the 114 A1w fish straight from the all-sites VCF, so the observed targets in METHODS.md
## section 3 are reproduced from the VCF rather than taken from the drift test's index_units.tsv.
## Every allele count must reproduce tail/index_units.tsv exactly.
##   bash code/11_Drift_vs_Selection/simulation/00_extract_index_genotypes.sh
set -euo pipefail
root="$(git rev-parse --show-toplevel)"
cd "${root}"
source "${root}/config.env"

work=${SIM_WORK:-${sweep_scratch}/peaks45_sim}
VCF="${sweep_scratch}/03_QC/${allsites_vcf}"
mkdir -p "${work}"
cd "${work}"

[ -f "${VCF}" ] || { echo "ERROR: all-sites VCF not found at ${VCF} (the 03_VCF_QC_Filter all-sites VCF from output_data/03_QC/; code/11_Drift_vs_Selection/sweep_detection/00_fetch_inputs.sb places it here)" >&2; exit 1; }

awk -F'\t' '$2=="BAM"{print $1}' "${root}/input_data/05_PopGen/${popsfile}" > bam.ids
awk -F'\t' '$2=="APA"||$2=="BEN"{print $1}' "${root}/input_data/05_PopGen/${popsfile}" > tp1.ids
cat bam.ids tp1.ids > bamtp1.ids
wc -l bam.ids tp1.ids

sx() { singularity exec --bind "${scratch_dir}" "${gwas_tools_image}" "$@"; }

sx bcftools view -r chr6:29503232,chr16:10916883,chr17:11692887 -S bamtp1.ids --force-samples "${VCF}" -Ov 2>/dev/null \
 | sx bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' - > idx_gt.raw
sx bcftools view -r chr6:29503232 -S bamtp1.ids --force-samples "${VCF}" -Ov 2>/dev/null \
 | sx bcftools query -l - > idx_gt.samples
wc -l idx_gt.raw idx_gt.samples
cut -f1,2,3,4 idx_gt.raw
exit 0
