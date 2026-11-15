## 04_submit_dinvestigate.sh
## August 2025
## JRG
## Submits one Dsuite Dinvestigate job per chromosome (run_dsuite_investigate.sb). Fig. S6 uses all 25 chromosomes
## (the commented loop); the active loop reruns only the six peak chromosomes.

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"config.env"

mkdir -p ${root}/output_data/slurm_logs/10_Introgression/

outdir=${root}/output_data/10_Introgression/
indir=${root}/input_data/10_Introgression/

mkdir -p ${outdir}

#for chr in {1..25}; do
 for chr in 13 16 17 24 6 8; do
	echo sbatch --job-name "DSUITE_chr${chr}" \
		   --output "${root}/output_data/slurm_logs/10_Introgression/dinvestigate_${chr}.log" \
		   --export=chr=${chr},root=${root} \
		   "${root}/code/10_Introgression/run_dsuite_investigate.sb"

	sbatch --job-name "DSUITE_chr${chr}" \
		   --output "${root}/output_data/slurm_logs/10_Introgression/dinvestigate_${chr}.log" \
		   --export=chr=${chr},root=${root} \
		   "${root}/code/10_Introgression/run_dsuite_investigate.sb"
done