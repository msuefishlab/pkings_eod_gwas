#!/bin/bash --login

## 04_Render_ADMIXTURE_report.sh
# Renders Admixture_Results.Rmd for one ADMIXTURE run (after 03_make_admixture_datafile.sh OUTNAME).
# Usage: bash code/13_Structure/04_Render_ADMIXTURE_report.sh OUTNAME

root="$(git rev-parse --show-toplevel)"

source $root/config.env

OUTNAME=$1

outdir=${root}/output_data/13_Structure/${OUTNAME}

mkdir -p ${outdir}

Rscript -e "rmarkdown::render('${root}/code/13_Structure/Admixture_Results.Rmd', output_dir = '${outdir}', output_file = '${OUTNAME}_admixture.report.html', params = list(outname = '${OUTNAME}'))"
