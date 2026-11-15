#!/bin/bash --login

OUTNAME=$1

MAF_THRESH=$2

# lrt (default; GEMMA -lmm 2, the main analysis) or wald (-lmm 1)
TEST=${3:-lrt}

root="$(git rev-parse --show-toplevel)"

source $root/config.env

outdir=${root}/output_data/06_Association/${OUTNAME}

mkdir -p ${outdir}

# Rscript ${scriptdir}/render_report.R INPUTDATA OUTPUTDATA REPORT.Rmd
singularity exec ${gwas_tools_image} Rscript ${root}/code/00_common/render_gwas_report.R ${OUTNAME} ${outdir}/${OUTNAME}_GWAS_Report ${root}/code/06_Association/gwas_report.Rmd ${MAF_THRESH} ${TEST}