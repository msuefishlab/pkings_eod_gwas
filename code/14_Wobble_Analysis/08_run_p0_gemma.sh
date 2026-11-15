#!/bin/bash

## 08_run_p0_gemma.sh
## Submits GEMMA Wald scans of each BAM_P0.pheno column (07_make_p0_phenos.R).
## Arg1 = array spec: "0" = real scans only; "1-100%25" = 100 permutations each.
## Optional Arg2 = column names to run (comma-separated; default all).
# e.g.  bash code/14_Wobble_Analysis/08_run_p0_gemma.sh 0
#       bash code/14_Wobble_Analysis/08_run_p0_gemma.sh "1-100%25" CLS_BIN,MANUAL_BIN

root="$(git rev-parse --show-toplevel)"
source ${root}/"config.env"

ARRAY=$1
ONLY=${2:-}
logdir=${root}/output_data/slurm_logs/14_Wobble_Analysis/BAM_P0
mkdir -p ${logdir}

tail -n +2 ${root}/input_data/14_Wobble_Analysis/BAM_P0/BAM_P0.pheno_columns.tsv | \
while IFS=$'\t' read -r COL NAME N NCASE NCTRL; do
  if [ -n "${ONLY}" ] && [[ ",${ONLY}," != *",${NAME},"* ]]; then continue; fi
  sbatch --job-name P0_${NAME} --array=${ARRAY} \
    --output ${logdir}/${NAME}_%a.slurm.log \
    --export=root=${root},COL=${COL},NAME=${NAME} \
    ${root}/code/14_Wobble_Analysis/p0_gemma.sb
done
