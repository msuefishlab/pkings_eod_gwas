## 04_run_dinvestigate.sh
## Dsuite Dinvestigate (-w 50,25, as in run_dsuite_investigate.sb) per chromosome for the six CURV3 trios in
## input_data/10_Introgression/test_trios_curv3.txt, CURV3 sets. Trios are written as in the main analysis (P1 = the
## lineage tested, P2 = its reference), so -f_dM = sharing between P1 and CURV3, the sign DsuiteInvestigateAnalysis.Rmd
## uses. Also reruns the main analysis's BP2-TP_OUT-BP_ANC trio on chr24 (test_trios_published_check.txt, SETS.txt) to
## check the windows against output_data/10_Introgression/BP2_TP_OUT_BP_ANC_localFstats_chr24_50_25.txt.
## Usage: bash 04_run_dinvestigate.sh [N_PARALLEL]
## Runs with Docker (Dsuite in the dsuite_arm image from 00_build_dsuite.sh).

#!/bin/bash

root="$(git rev-parse --show-toplevel)"

npar=${1:-7}
rel=output_data/10_Introgression/peak6_curv3
scratch=${root}/${rel}/_scratch
mkdir -p ${root}/${rel}/dinvestigate_curv3 ${root}/${rel}/dinvestigate_published_check ${root}/${rel}/logs

run_chr() {
	set_name=$1; chr=$2; root=$3; rel=$4
	if [ ${set_name} == curv3 ]; then sets=SETS_curv3.txt; trios=test_trios_curv3.txt; outdir=dinvestigate_curv3
	else sets=SETS.txt; trios=test_trios_published_check.txt; outdir=dinvestigate_published_check; fi
	docker run --rm -v ${root}:/project_root -w /project_root/${rel}/${outdir} dsuite_arm \
		/opt/Dsuite/Build/Dsuite Dinvestigate -w 50,25 -n chr${chr} \
		/project_root/${rel}/_scratch/chr_vcfs/chr${chr}.snps.gt.vcf.gz \
		/project_root/input_data/10_Introgression/${sets} \
		/project_root/input_data/10_Introgression/${trios} < /dev/null \
		> ${root}/${rel}/logs/dinvestigate_${set_name}_chr${chr}.log 2>&1
	if grep -q "^chr" ${root}/${rel}/${outdir}/*_localFstats_chr${chr}_50_25.txt 2>/dev/null; then
		echo "${set_name} chr${chr} finished $(date +%T)"; else echo "${set_name} chr${chr} FAILED"; fi
}
export -f run_chr

# Stop unless all 25 chromosome VCFs from 02_make_chr_vcfs.sh are there
for chr in {1..25}; do
	[ -s ${scratch}/chr_vcfs/chr${chr}.snps.gt.vcf.gz.tbi ] || { echo "missing chr${chr} VCF: run 02_make_chr_vcfs.sh"; exit 1; }
done

# Largest chromosomes (most SNPs) first; the published-trio check last
chrs=$(for chr in {1..25}; do echo "$(bcftools index -n ${scratch}/chr_vcfs/chr${chr}.snps.gt.vcf.gz) ${chr}"; done | sort -k1,1nr | cut -d' ' -f2)
{ for chr in ${chrs}; do echo "curv3 ${chr}"; done; echo "published 24"; } |
	xargs -P ${npar} -L 1 bash -c 'run_chr $0 $1 '"${root} ${rel}"
