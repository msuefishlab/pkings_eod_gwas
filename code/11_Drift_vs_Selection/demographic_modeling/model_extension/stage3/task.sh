#!/bin/bash
## task.sh '<out_check> <tab> <kind> <args...>' -- one idempotent Stage 3 task line; kind: sim (sim_task.py), post
## (genome_post.py), rec (records_full.py), obs (obs_rows.py). Skips the task when <out_check> exists and is non-empty.
set -eu
IFS=$'\t' read -r chk kind rest <<< "$1"
[ -s "$chk" ] && exit 0
s=$(date +%s)
case "$kind" in
  sim)  $PY $X/sim_task.py $rest ;;
  post) $PY $X/genome_post.py $rest ;;
  rec)  $PY $X/records_full.py $rest ;;
  obs)  $PY $X/obs_rows.py $rest ;;
  *) echo "bad kind $kind"; exit 2 ;;
esac
echo "WALL $kind $(( $(date +%s) - s ))s $chk"
