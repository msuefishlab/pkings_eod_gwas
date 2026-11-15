#!/bin/bash

## 05_collect_min_values.sh
## Feb 2023
## JRG
# This script collects minimum p-values across permuted gemma gwas runs, which would normally be collected by gemma-wrapper when not in slurm mode.
# e.g. . code/06_Association/05_collect_min_values.sh PKINGS_ALL_WOB_EXCLUDED

root="$(git rev-parse --show-toplevel)"


OUTNAME=$1

cd $root

source $root/config.env

indir=${root}/output_data/06_Association/
outdir=${root}/output_data/06_Association/$OUTNAME/${OUTNAME}_permution/

rm -rf ${outdir}/min.txt
touch ${outdir}/min.txt
# the p-value column is found by name (p_wald for -lmm 1, p_lrt for -lmm 2); nan/-nan rows are skipped
for f in ${outdir}/*.assoc.txt;
	do
		if grep -e "phenotypes-" ${f%%.assoc.txt}.log.txt;
			then echo "Processing $f file...";
				awk 'NR==1 { for (i = 1; i <= NF; i++) if ($i == "p_wald" || $i == "p_lrt") pcol = i; next }
				     $pcol ~ /^[0-9.eE+-]+$/ && (!found || $pcol + 0 < min) { line = $0; min = $pcol + 0; found = 1 }
				     END { print line }' $f >> ${outdir}/min.txt;
		fi;
	done
echo "$(wc -l < ${outdir}/min.txt) permutations collected into ${outdir}/min.txt"