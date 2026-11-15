#!/bin/bash
## run_tasks_node.sh <tasks.txt> [<post.txt>] -- one node: unpack the records cache (node-local $TMPDIR when it has
## > 150 GB free, else $O/cache), run every task line through task.sh with xargs (idempotent: finished outputs are
## skipped), then the post lines. Logs to <tasks>.log.
set -euo pipefail
source "$(dirname "$0")/env.sh"
NP=$(( ${SLURM_CPUS_PER_TASK:-64} - 2 ))
T=$1; PT=${2:-}
free=$(df -BG --output=avail ${TMPDIR:-/tmp} | tail -1 | tr -dc 0-9)
if [ "${free:-0}" -gt 150 ]; then CACHE=${TMPDIR}/s3cache; else CACHE=$O/cache; fi
echo "host $(hostname) cpus ${SLURM_CPUS_PER_TASK:-} workers $NP cache $CACHE start $(date)"
$PY $X/prep_cache.py $O/records $S3_MASK $CACHE 64
export S3_CACHE=$CACHE
echo "tasks $(grep -c . $T) $(date)"
grep . $T | xargs -P $NP -I{} bash $X/task.sh '{}' > $T.log 2>&1 || echo "WARN some tasks failed"
if [ -n "$PT" ]; then
  echo "post $(grep -c . $PT) $(date)"
  grep . $PT | xargs -P $NP -I{} bash $X/task.sh '{}' > $PT.log 2>&1 || echo "WARN some post tasks failed"
fi
n_left=0
while IFS=$'\t' read -r chk rest; do [ -s "$chk" ] || n_left=$((n_left+1)); done < <(grep . $T; [ -n "$PT" ] && grep . $PT)
echo "unfinished outputs: $n_left $(date)"
exit 0
