#!/bin/bash
## orch_genomes.sh <code_dir> <tag> <i_ratio> <model:ratio[,model:ratio...]> <n_genomes> -- whole M-sum replicate
## genomes g = 1..n for every listed model at its r/mu (sim_index ratio index i_ratio), each through ext_aggregate
## (genome_post.py); returns one tarball of the per-genome panel.tsv / ld_curves.tsv / chromld.npz.
set -uo pipefail
source $1/env.sh
TAG=$2; IR=$3; LIST=$4; NG=$5
JOBDIR=$PWD
cd $O; mkdir -p genomes logs
TF=$O/genomes/tasks_$TAG.txt; PF=$O/genomes/post_$TAG.txt
if [ ! -s $TF ]; then
  : > $TF.tmp; : > $PF.tmp
  for mr in ${LIST//,/ }; do $PY $X/gen_tasks.py genomes ${mr%%:*} $IR ${mr##*:} 1 $NG $TF.tmp $PF.tmp; done
  mv $PF.tmp $PF; mv $TF.tmp $TF
fi
bash $X/submit_node.sh s3_$TAG 23:00:00 $X/run_tasks_node.sh $TF $PF; echo "node exit $?"
cd $O/genomes
files=""
for mr in ${LIST//,/ }; do files="$files $(ls -d ${mr%%:*}/r$IR/g*/panel.tsv ${mr%%:*}/r$IR/g*/ld_curves.tsv ${mr%%:*}/r$IR/g*/chromld.npz 2>/dev/null | tr '\n' ' ')"; done
tar czf $O/genomes/summary_$TAG.tgz $files
cd $JOBDIR
cp $O/genomes/summary_$TAG.tgz ./
grep -E 'Error|Traceback' $TF.log $PF.log 2>/dev/null | head -20 > task_errors.txt || true
grep 'WALL sim' $TF.log | awk '{print $3}' | sort -n | awk '{a[NR]=$1} END{print "n", NR, "median", a[int(NR/2)+1], "max", a[NR]}' > walls.txt 2>/dev/null || true
cp $(ls -t $O/logs/s3_${TAG}_*.out | head -1) ./node.out 2>/dev/null || true
ls -la
exit 0
