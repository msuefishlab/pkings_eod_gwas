#!/bin/bash

## 09_reduce_p0_gemma.sh
## On the cluster, after 08_run_p0_gemma.sh: shrinks each BAM_P0 scan to what the
## summary (09_collect_p0_gemma.R) needs, so the ~500 MB assoc files need not be copied.
##   <NAME>.top.tsv       SNPs with p_wald < 1e-3 (header kept)
##   <NAME>.perm_min.tsv  minimum-p line of every finished permutation
# e.g.  sbatch --wrap "bash code/14_Wobble_Analysis/09_reduce_p0_gemma.sh"
#       then copy output_data/14_Wobble_Analysis/BAM_P0/*/*.{top,perm_min}.tsv to the machine that runs 09_collect_p0_gemma.R

root="$(git rev-parse --show-toplevel)"
base=${root}/output_data/14_Wobble_Analysis/BAM_P0

for d in ${base}/*/; do
  NAME=$(basename ${d})
  [ -f ${d}/${NAME}.lmm1.assoc.txt ] || continue
  # column 12 = p_wald
  awk 'NR==1 || $12 < 1e-3' ${d}/${NAME}.lmm1.assoc.txt > ${d}/${NAME}.top.tsv
  first=$(ls ${d}/perm/perm_*.min.txt 2>/dev/null | head -1)
  if [ -n "${first}" ]; then
    head -1 ${first} | awk '{print "perm\t"$0}' > ${d}/${NAME}.perm_min.tsv
    for f in ${d}/perm/perm_*.min.txt; do
      k=$(basename ${f} .min.txt); k=${k#perm_}
      tail -n +2 ${f} | awk -v k=${k} '{print k"\t"$0}' >> ${d}/${NAME}.perm_min.tsv
    done
  fi
  echo "${NAME}: $(($(wc -l < ${d}/${NAME}.top.tsv) - 1)) SNPs p<1e-3; $(ls ${d}/perm/perm_*.min.txt 2>/dev/null | wc -l) perms"
done
