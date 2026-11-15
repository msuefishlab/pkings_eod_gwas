#!/bin/bash
## orch_step3.sh <code_dir> <models,comma> <tag> -- Stage 3 step 3 simulations: reference families (cal, ins) and the
## held-out scoring family for every listed model x r/mu grid value of simulation/23e; pooled into fam_summary_<tag>.npz.
set -uo pipefail
source $1/env.sh
MODELS=$2; TAG=$3
cd $O; mkdir -p fam logs
TF=$O/fam/tasks_$TAG.txt
[ -s $TF ] || $PY $X/gen_tasks.py families $MODELS $TF
bash $X/submit_node.sh s3_step3_$TAG 23:00:00 $X/run_tasks_node.sh $TF; echo "node exit $?"
$PY $X/family_summary.py $O/fam $MODELS $O/fam/fam_summary_$TAG.npz > logs/fam_summary_$TAG.log 2>&1
cat logs/fam_summary_$TAG.log
cd - >/dev/null
cp $O/fam/fam_summary_$TAG.npz $O/fam/draws.tsv ./ 2>/dev/null || true
grep -E 'Error|Traceback' $TF.log | head -20 > task_errors.txt 2>/dev/null || true
grep -c WALL $TF.log > n_done.txt 2>/dev/null || true
cp $(ls -t $O/logs/s3_step3_${TAG}_*.out | head -1) ./node.out 2>/dev/null || true
ls -la
exit 0
