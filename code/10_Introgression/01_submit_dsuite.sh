## 01_submit_dsuite.sh
## August 2025
## JRG
## Submits one Dsuite Dtrios job per chromosome (run_dsuite.sb). Combine the 25 outputs with DtriosCombine (README).

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"config.env"

mkdir -p ${root}/output_data/slurm_logs/10_Introgression/

outdir=${root}/output_data/10_Introgression/
indir=${root}/input_data/10_Introgression/

mkdir -p ${outdir}

for chr in {1..25}; do
	echo sbatch --job-name "DSUITE_chr${chr}" \
		   --output "${root}/output_data/slurm_logs/10_Introgression/dsuite_${chr}.log" \
		   --export=chr=${chr},root=${root} \
		   "${root}/code/10_Introgression/run_dsuite.sb"

	sbatch --job-name "DSUITE_chr${chr}" \
		   --output "${root}/output_data/slurm_logs/10_Introgression/dsuite_${chr}.log" \
		   --export=chr=${chr},root=${root} \
		   "${root}/code/10_Introgression/run_dsuite.sb"
done