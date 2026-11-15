#!/bin/bash --login

OUTNAME=$1

MAF_THRESH=$2

root="$(git rev-parse --show-toplevel)"

source $root/config.env

outdir=${root}/output_data/08_Peak_Analysis/${OUTNAME}

mkdir -p ${outdir}

# Rscript ${scriptdir}/render_report.R INPUTDATA OUTPUTDATA REPORT.Rmd
singularity exec ${gwas_tools_image} Rscript ${root}/code/00_common/render_rmd_report.R ${OUTNAME} ${outdir}/${OUTNAME}_PVE_Report ${root}/code/08_Peak_Analysis/PVE_new.Rmd ${MAF_THRESH}