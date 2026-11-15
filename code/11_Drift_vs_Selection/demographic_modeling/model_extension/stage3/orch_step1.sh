#!/bin/bash
## orch_step1.sh -- controlling job, Stage 3 step 1. Writes sim_index.tsv (before any simulation), then runs
## step1_node.sh on one node (submit_node.sh) and copies the small result tables into the job workdir.
set -uo pipefail
source $1/env.sh
mkdir -p $O/logs
cd $O
if [ ! -s sim_index.tsv ]; then
  $PY $X/write_sim_index.py $S3_PARAMS $FDA/fits/best_fits.tsv sim_index.tsv
fi
bash $X/submit_node.sh s3_step1 12:00:00 $X/step1_node.sh; echo "node exit $?"
cd - >/dev/null
for f in sim_index.tsv records_checks.tsv records_summary.txt dropout_d.tsv pilot_repro.tsv pilot_repro_detail.tsv RECORDS_MISMATCH STEP1.DONE logs/provenance_step1.txt logs/pilot_compare.log logs/dropout.log; do
  [ -e $O/$f ] && cp $O/$f ./$(basename $f) || true
done
ls $O/logs/s3_step1_*.out >/dev/null 2>&1 && cp $(ls -t $O/logs/s3_step1_*.out | head -1) ./node_step1.out || true
cp $O/test_genome/M2/r99/g1/panel.tsv ./test_genome_panel.tsv 2>/dev/null || true
cp $O/test_genome/M2/r99/g1/ld_curves.tsv ./test_genome_ld_curves.tsv 2>/dev/null || true
grep -E 'WALL|cpu=|Error|error' $O/logs/test_genome.log > ./test_genome_log.txt 2>/dev/null || true
grep -E 'WALL' $O/logs/pilot_repro.log | awk '{print $3}' | sort -n | tail -3 > ./pilot_task_walls.txt 2>/dev/null || true
grep -E 'Error|error|Traceback' $O/logs/records.log | head -20 > ./records_errors.txt 2>/dev/null || true
tail -5 $O/logs/records.log > ./records_log_tail.txt 2>/dev/null || true
ls -la
exit 0
