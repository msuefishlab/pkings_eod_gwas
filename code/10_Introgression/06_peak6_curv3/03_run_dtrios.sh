## 03_run_dtrios.sh
## Dsuite Dtrios with the CURV3 sets (input_data/10_Introgression/SETS_curv3.txt, pops_curv3.tree):
##   1. the extended Peak 6 interval, cut from the VCF exactly as in 01_validate_peak6.sh (all samples);
##   2. genome-wide, one Dtrios per chromosome on the 02_make_chr_vcfs.sh files, then DtriosCombine (as in the main
##      analysis: code/10_Introgression/run_dsuite.sb, whole_genome_combined_*).
## It also reruns the main sets (SETS.txt, pops.tree) genome-wide on the same files, to check the genome-wide baseline
## against the main analysis's whole_genome_combined_tree.txt. CURV3 jobs are queued first.
## Usage: bash 03_run_dtrios.sh [VCF] [N_PARALLEL]
## Runs with Docker (Dsuite in the dsuite_arm image from 00_build_dsuite.sh).

#!/bin/bash

root="$(git rev-parse --show-toplevel)"

vcf=${1:-${root}/output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz}
npar=${2:-7}
rel=output_data/10_Introgression/peak6_curv3
scratch=${root}/${rel}/_scratch
mkdir -p ${root}/${rel}/peak_curv3 ${root}/${rel}/genome_curv3 ${root}/${rel}/genome_published_sets ${root}/${rel}/logs

dsuite() {
	docker run --rm -v ${root}:/project_root -w /project_root dsuite_arm /opt/Dsuite/Build/Dsuite "$@" < /dev/null
}

## 1. Peak 6, CURV3 sets
read -r chr start end peak _ < <(awk '$4 == 6' ${root}/input_data/10_Introgression/peaks.bed)
if [ ! -s ${root}/${rel}/peak_curv3/peak_${peak}_tree.txt ]; then
	bcftools view -r chr${chr}:${start}-${end} -v snps -m2 -M2 ${vcf} -Oz -o ${scratch}/peak_${peak}_curv3.vcf.gz
	dsuite Dtrios -o /project_root/${rel}/peak_curv3/peak_${peak} \
		-t /project_root/input_data/10_Introgression/pops_curv3.tree \
		/project_root/${rel}/_scratch/peak_${peak}_curv3.vcf.gz \
		/project_root/input_data/10_Introgression/SETS_curv3.txt > ${root}/${rel}/logs/dtrios_peak_${peak}_curv3.log 2>&1
	rm -f ${scratch}/peak_${peak}_curv3.vcf.gz
fi

## 2. Genome-wide, per chromosome (outputs genome_curv3/ and genome_published_sets/; finished chromosomes are skipped)
run_chr() {
	set_name=$1; chr=$2; root=$3; rel=$4
	if [ ${set_name} == curv3 ]; then sets=SETS_curv3.txt; tree=pops_curv3.tree; else sets=SETS.txt; tree=pops.tree; fi
	out=${rel}/genome_${set_name}/chr${chr}
	if [ -s ${root}/${out}_tree.txt ]; then echo "${set_name} chr${chr} done already"; return 0; fi
	docker run --rm -v ${root}:/project_root -w /project_root dsuite_arm /opt/Dsuite/Build/Dsuite Dtrios \
		-o /project_root/${out} -t /project_root/input_data/10_Introgression/${tree} \
		/project_root/${rel}/_scratch/chr_vcfs/chr${chr}.snps.gt.vcf.gz \
		/project_root/input_data/10_Introgression/${sets} < /dev/null \
		> ${root}/${rel}/logs/dtrios_${set_name}_chr${chr}.log 2>&1
	if [ -s ${root}/${out}_tree.txt ]; then echo "${set_name} chr${chr} finished $(date +%T)"; else echo "${set_name} chr${chr} FAILED"; fi
}
export -f run_chr

# Stop unless all 25 chromosome VCFs from 02_make_chr_vcfs.sh are there
for chr in {1..25}; do
	[ -s ${scratch}/chr_vcfs/chr${chr}.snps.gt.vcf.gz.tbi ] || { echo "missing chr${chr} VCF: run 02_make_chr_vcfs.sh"; exit 1; }
done

# Largest chromosomes (most SNPs) first
chrs=$(for chr in {1..25}; do echo "$(bcftools index -n ${scratch}/chr_vcfs/chr${chr}.snps.gt.vcf.gz) ${chr}"; done | sort -k1,1nr | cut -d' ' -f2)
{ for chr in ${chrs}; do echo "curv3 ${chr}"; done; for chr in ${chrs}; do echo "published_sets ${chr}"; done; } |
	xargs -P ${npar} -L 1 bash -c 'run_chr $0 $1 '"${root} ${rel}"

## 3. Combine the chromosomes (DtriosCombine), as for the main analysis's whole_genome_combined_*
for set_name in curv3 published_sets; do
	if [ ${set_name} == curv3 ]; then tree=pops_curv3.tree; else tree=pops.tree; fi
	dsuite DtriosCombine -o /project_root/${rel}/genome_${set_name}/whole_genome \
		-t /project_root/input_data/10_Introgression/${tree} \
		$(for chr in {1..25}; do echo /project_root/${rel}/genome_${set_name}/chr${chr}; done) \
		> ${root}/${rel}/logs/dtrios_combine_${set_name}.log 2>&1
done
