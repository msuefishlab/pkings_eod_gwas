#!/bin/bash

## 02_create_sample_batches.sh
## Feb 2023
## JRG
# This script splits the samples of the SNP VCF into batches of 25 (input_data/07_Phasing/sample_names_bNN) and is
# run interactively from the command line.

root="$(git rev-parse --show-toplevel)"
source ${root}/"config.env"

outdir=${root}/input_data/07_Phasing/
mkdir -p ${outdir}

module purge; module load BCFtools

bcftools query -l output_data/03_QC/${snps_only_for_phasing} | split -l 25 -d - input_data/07_Phasing/sample_names_b

