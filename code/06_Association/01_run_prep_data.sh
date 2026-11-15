#!/bin/bash

## 01_run_prep_data.sh
## January 2024
## JRG
# This script prepares a VCF files for a set of individuals specified in arg1 with phenotypes in arg2 for relationship matrix construction and association analysis.
# e.g. bash code/06_Association/01_run_prep_data.sh PKINGS_ALL_WOB_EXCLUDED

root="$(git rev-parse --show-toplevel)"
source ${root}/"config.env"

mkdir -p ${root}/output_data/slurm_logs/06_Assocation/

OUTNAME=$1

echo sbatch --job-name $OUTNAME"_PREP_DATA" --output ${root}"/output_data/slurm_logs/06_Assocation/ASSOC_"${OUTNAME}"_prep_data.log" --export=root=${root},OUTNAME=${OUTNAME} ${root}/code/06_Association/prep_data.sb
sbatch --job-name $OUTNAME"_PREP_DATA" --output ${root}"/output_data/slurm_logs/06_Assocation/ASSOC_"${OUTNAME}"_prep_data.log" --export=root=${root},OUTNAME=${OUTNAME} ${root}/code/06_Association/prep_data.sb
