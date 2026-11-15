#!/bin/bash

## 04_squish_outputs.sh
## Feb 2023
## JRG
# This script merges each sample's WhatsHap outputs for the unplaced scaffolds into one final_chunk VCF
# e.g. bash code/07_Phasing/04_squish_outputs.sh input_data/07_Phasing/sample_names_b00

root="$(git rev-parse --show-toplevel)"
source ${root}/"config.env"

outdir=${root}/input_data/07_Phasing/
mkdir -p ${outdir}

mkdir -p ${root}/output_data/slurm_logs/07_Phasing/

while read SAMPLE
do
  echo "$SAMPLE"
  mkdir -p ${outdir}/$SAMPLE

  echo sbatch --job-name ${SAMPLE}_WHATSHAP --output ${root}/output_data/slurm_logs/07_Phasing/${SAMPLE}_SQUISH_OUTPUTS.slurm.log --export=root=${root} ${root}/code/07_Phasing/squish_outputs.sb $SAMPLE
  sbatch --job-name ${SAMPLE}_WHATSHAP --output ${root}/output_data/slurm_logs/07_Phasing/${SAMPLE}_SQUISH_OUTPUTS.slurm.log --export=root=${root} ${root}/code/07_Phasing/squish_outputs.sb $SAMPLE

done < "$1"