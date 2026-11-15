## 03_submit_dsuite_peaks.sh
## August 2025
## JRG
## Submits Dsuite Dtrios for each GWAS peak window.
##
## Windows are the extended peaks in input_data/10_Introgression/peaks.bed (-log10 p > 3 merged within 75 kb), used
## as-is.
##
## Input is the unphased SNP VCF the genome-wide run used (03_QC/...FILTERED_BIALLELIC_SNPS_ONLY.renamed.vcf.gz).
## Pass another path as $1 if needed; it must be bgzipped and indexed.

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"config.env"

vcf=${1:-${root}/output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_SNPS_ONLY.renamed.vcf.gz}

mkdir -p ${root}/output_data/slurm_logs/10_Introgression/
mkdir -p ${root}/output_data/10_Introgression/peaks_extended/

while read -r chr start end peak _; do

	sbatch --job-name "DSUITE_peak${peak}" \
		   --output "${root}/output_data/slurm_logs/10_Introgression/dsuite_peaks_peak_${peak}.log" \
		   --export=chr=${chr},root=${root},start=${start},end=${end},peak=${peak},vcf=${vcf} \
		   "${root}/code/10_Introgression/run_dsuite_peaks.sb"

done < ${root}/input_data/10_Introgression/peaks.bed
