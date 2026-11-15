## 02_make_chr_vcfs.sh
## Cuts per-chromosome biallelic SNP VCFs (chr1-chr25) from the SNP + indel VCF for the genome-wide Dtrios and
## Dinvestigate runs: `bcftools view -v snps -m2 -M2`, the 209 samples in SETS.txt (the CURV3 sets use the same
## samples), GT only. Dsuite reads only GT and ignores samples outside the SETS file, and dropping samples, INFO and
## other FORMAT fields keeps every record, so the line count (Dsuite's jackknife blocks) is unchanged. The files are
## a few GB in all and are deleted by 06_cleanup.sh.
## Usage: bash 02_make_chr_vcfs.sh [VCF] [N_PARALLEL]
## Runs with bcftools and tabix.

#!/bin/bash

root="$(git rev-parse --show-toplevel)"

vcf=${1:-${root}/output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz}
npar=${2:-8}
scratch=${root}/output_data/10_Introgression/peak6_curv3/_scratch
mkdir -p ${scratch}/chr_vcfs

cut -f1 ${root}/input_data/10_Introgression/SETS.txt > ${scratch}/sets_samples.txt

make_one() {
	chr=$1; vcf=$2; scratch=$3
	out=${scratch}/chr_vcfs/chr${chr}.snps.gt.vcf.gz
	if [ -s ${out}.tbi ]; then echo "chr${chr} done already"; return 0; fi
	bcftools view -r chr${chr} -v snps -m2 -M2 -S ${scratch}/sets_samples.txt -Ou ${vcf} 2>/dev/null |
		bcftools annotate -x INFO,^FORMAT/GT -Oz -o ${out}.tmp &&
		mv ${out}.tmp ${out} && tabix -p vcf ${out} &&
		echo "chr${chr} $(bcftools index -n ${out}) SNPs"
}
export -f make_one

# All 25 chromosomes, N_PARALLEL at a time
# (BSD xargs -I caps the command at 255 bytes, so the arguments are appended with -L 1 instead)
for chr in {1..25}; do echo ${chr}; done |
	xargs -P ${npar} -L 1 bash -c 'make_one $0 '"${vcf} ${scratch}"
