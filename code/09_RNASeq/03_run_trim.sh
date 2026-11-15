#!/bin/bash

## Run by bash code/09_RNASeq/03_run_trim.sh input_data/09_RNASeq/rnaseq_metadata.txt
## Trimmomatic on output_data/09_RNASeq/raw_reads/<SpecNo>_R{1,2}_combined.fastq.gz, one job per specimen (run_trim.sb).

root="$(git rev-parse --show-toplevel)"
source "$root/config.env"

# Define the input file
input_file=$1

readsdir="${root}/output_data/09_RNASeq/raw_reads"  # Define the output director
outdir="${root}/output_data/09_RNASeq/trimmed_reads"  # Define the output director
slurmoutdir=${root}/output_data/slurm_outputs/09_RNASeq

mkdir -p $outdir
mkdir -p $slurmoutdir

{
  read -r header # Skip the header
  declare -A seen_specnos # Declare an associative array to track unique specnos

  while IFS=$'\t' read -r specno population r1 r2 phenotype; do
    if [[ -z "${seen_specnos[$specno]}" ]]; then
      # If specno is not yet seen
      seen_specnos[$specno]=1 # Mark specno as seen
      mkdir -p "${outdir}/${specno}"
      echo sbatch --output="${slurmoutdir}/trim_${specno}.out" --export=READ1="${readsdir}/${specno}_R1_combined.fastq.gz",READ2="${readsdir}/${specno}_R2_combined.fastq.gz",OUTNAME="${outdir}/${specno}",SAMPLE="${specno}" "${root}/code/09_RNASeq/run_trim.sb"
      sbatch --output="${slurmoutdir}/trim_${specno}.out" --export=READ1="${readsdir}/${specno}_R1_combined.fastq.gz",READ2="${readsdir}/${specno}_R2_combined.fastq.gz",OUTNAME="${outdir}/${specno}",SAMPLE="${specno}" "${root}/code/09_RNASeq/run_trim.sb"
    fi
  done
} < "$input_file"

