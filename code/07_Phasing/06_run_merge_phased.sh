#!/bin/bash

## 06_run_merge_phased.sh
## Feb 2023
## JRG
# This script concatenates the phased regions into one VCF
# e.g. bash code/07_Phasing/06_run_merge_phased.sh

root="$(git rev-parse --show-toplevel)"
source ${root}/"config.env"

outdir=${root}/input_data/07_Phasing/
mkdir -p ${outdir}

mkdir -p ${root}/output_data/slurm_logs/07_Phasing/

echo sbatch --job-name "MERGE_PHASED" --output ${root}/output_data/slurm_logs/07_Phasing/MERGE_PHASED.slurm.log --export=root=${root} ${root}/code/07_Phasing/merge_phased.sb
sbatch --job-name "MERGE_PHASED" --output ${root}/output_data/slurm_logs/07_Phasing/MERGE_PHASED.slurm.log --export=root=${root} ${root}/code/07_Phasing/merge_phased.sb