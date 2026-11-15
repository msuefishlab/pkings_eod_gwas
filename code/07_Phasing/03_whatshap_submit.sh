#!/bin/bash

## 03_whatshap_submit.sh
## Feb 2023
## JRG
# This script runs WhatsHap read-backed phasing for every sample in a batch file, one array job per region list
# e.g. bash code/07_Phasing/03_whatshap_submit.sh input_data/07_Phasing/sample_names_b00

root="$(git rev-parse --show-toplevel)"
source ${root}/"config.env"

outdir=${root}/input_data/07_Phasing/
mkdir -p ${outdir}

mkdir -p ${root}/output_data/slurm_logs/07_Phasing/

while read SAMPLE
do
  echo "$SAMPLE"
  mkdir -p ${outdir}/$SAMPLE

  echo sbatch --job-name ${SAMPLE}_WHATSHAP --array=0-25 --output ${root}/output_data/slurm_logs/07_Phasing/${SAMPLE}_WHATSHAP.slurm_%a.log --export=root=${root} ${root}/code/07_Phasing/submit_whatshap.sb $SAMPLE
  sbatch --job-name ${SAMPLE}_WHATSHAP --array=0-25 --output ${root}/output_data/slurm_logs/07_Phasing/${SAMPLE}_WHATSHAP.slurm_%a.log --export=root=${root} ${root}/code/07_Phasing/submit_whatshap.sb $SAMPLE

done < "$1"