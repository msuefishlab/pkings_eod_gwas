## 00_copy_ar_bams_from_terra.sh
## September 2024
## JRG
## Copies the analysis-ready BAMs written by the 01_Terra workflow from the workspace bucket to ${scratch_dir}.
## input_data/07_Phasing/AR_BAMLIST_2024.txt is not shipped: write it yourself, one bucket URL per line for the BAMs
## in your ${bucket} (the analysis-ready BAM column of your Terra sample table). Skip this step if your BAMs are
## already on the cluster.

#!/bin/bash

root="$(git rev-parse --show-toplevel)"
source ${root}/"config.env"

mkdir -p ${root}/output_data/slurm_logs/07_Phasing/

readarray -t readfiles < ${root}/input_data/07_Phasing/AR_BAMLIST_2024.txt 

nbams=${#readfiles[@]}

echo sbatch --job-name "COPY_AR_BAMS" --output ${root}"/output_data/slurm_logs/07_Phasing/PHASING_COPY_AR_BAMS_%a.log" -a 0-$nbams --export=root=${root} ${root}/code/07_Phasing/copy_ar_bams.sb
sbatch --job-name "COPY_AR_BAMS" --output ${root}"/output_data/slurm_logs/07_Phasing/PHASING_COPY_AR_BAMS_%a.log" -a 0-$nbams --export=root=${root} ${root}/code/07_Phasing/copy_ar_bams.sb
