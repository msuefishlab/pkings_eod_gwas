## 01_validate_peak6.sh
## Validation gate: rerun the main trio set (input_data/10_Introgression/SETS.txt and pops.tree) on the extended Peak 6
## interval, then compare with the Peak 6 rows of Supplementary Data S7 (validate_vs_dataS7.R). No CURV3 result is
## read unless this reproduces Data S7 (BP2-BP_ANC D = 0.197, 445 informative sites).
##
## The main analysis cut the region from the unphased biallelic SNP VCF (03_QC/...FILTERED_BIALLELIC_SNPS_ONLY.renamed.vcf.gz),
## all samples. This script reads ...FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz, the same SNP shards concatenated with
## the indel shards (code/03_VCF_QC_Filter/05_combine_filtered.sh), so `-v snps -m2 -M2` recovers the SNP set. All
## samples are kept (Dsuite ignores samples not in SETS.txt). Pass another VCF path as $1 if needed.
## Runs with bcftools and Docker (Dsuite in the dsuite_arm image from 00_build_dsuite.sh).

#!/bin/bash

root="$(git rev-parse --show-toplevel)"

vcf=${1:-${root}/output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz}
scratch=${root}/output_data/10_Introgression/peak6_curv3/_scratch
outdir=${root}/output_data/10_Introgression/peak6_curv3/validation_published_sets
mkdir -p ${scratch} ${outdir}

# Extended Peak 6, used exactly as 03_submit_dsuite_peaks.sh used peaks.bed
read -r chr start end peak _ < <(awk '$4 == 6' ${root}/input_data/10_Introgression/peaks.bed)

bcftools view -r chr${chr}:${start}-${end} -v snps -m2 -M2 ${vcf} -Oz -o ${scratch}/peak_${peak}_published.vcf.gz

docker run --rm -v ${root}:/project_root -w /project_root dsuite_arm \
	/opt/Dsuite/Build/Dsuite Dtrios \
	-o /project_root/output_data/10_Introgression/peak6_curv3/validation_published_sets/peak_${peak} \
	-t /project_root/input_data/10_Introgression/pops.tree \
	/project_root/output_data/10_Introgression/peak6_curv3/_scratch/peak_${peak}_published.vcf.gz \
	/project_root/input_data/10_Introgression/SETS.txt < /dev/null

rm -f ${scratch}/peak_${peak}_published.vcf.gz

Rscript ${root}/code/10_Introgression/06_peak6_curv3/validate_vs_dataS7.R
