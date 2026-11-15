#!/bin/bash
## alt1_node.sh <models,comma> -- SFS unit-check simulations (sim_task.py sfs) for the listed models, i_ratio 40,
## physical r/mu 1.0, the 100 in-sample-family blocks x g = 1, 2 (200 block-sims per model, 4 tasks per model).
set -euo pipefail
source "$(dirname "$0")/env.sh"
NP=$(( ${SLURM_CPUS_PER_TASK:-64} - 2 ))
free=$(df -BG --output=avail ${TMPDIR:-/tmp} | tail -1 | tr -dc 0-9)
if [ "${free:-0}" -gt 150 ]; then CACHE=${TMPDIR}/s3cache; else CACHE=$O/cache; fi
$PY $X/prep_cache.py $O/records $S3_MASK $CACHE 64
export S3_CACHE=$CACHE
mkdir -p $O/alt/sfs/lists
awk -F'\t' 'NR>1 && $1=="ins"{print $4}' $O/fam/draws.tsv > $O/alt/sfs/lists/blocks.txt
: > $O/alt/sfs/tasks.txt
for m in ${1//,/ }; do for p in 0 1 2 3; do
  awk -v p=$p 'NR%4==p {print 1"\t"$1; print 2"\t"$1}' $O/alt/sfs/lists/blocks.txt > $O/alt/sfs/lists/${m}_p$p.tsv
  printf '%s/alt/sfs/%s_p%s.npz\tsim\tsfs %s 40 1.0 %s/alt/sfs/lists/%s_p%s.tsv %s/alt/sfs/%s_p%s.npz\n' $O $m $p $m $O $m $p $O $m $p >> $O/alt/sfs/tasks.txt
done; done
grep . $O/alt/sfs/tasks.txt | xargs -P $NP -I{} bash $X/task.sh '{}' > $O/alt/sfs/tasks.log 2>&1 || echo "WARN sfs tasks failed"
exit 0
