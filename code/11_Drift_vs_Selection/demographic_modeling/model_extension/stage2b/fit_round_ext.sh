#!/bin/bash
# One Stage 2b fitting round on one node: every line of $1 is a full fit_ext.py argument list (07 args + --bounds
# --theta_ref). Same structure as the drift test's fitting rounds. Idempotent: a task whose seed file already holds header +
# all starts is skipped.
set -uo pipefail
TASKS=$1; FDABS=$2; CODE=$3; IMG=$4
cd "$(dirname "$FDABS")"
R=${R:-${root:-$(git rev-parse --show-toplevel)}}
[ -n "${scratch_dir:-}" ] || { cd "$R" && source ./config.env && cd - >/dev/null; }
export OMP_NUM_THREADS=1 FDABS CODE IMG R scratch_dir
one() {
  args="$1"
  pair=$(sed -E 's/.*--pair ([^ ]+).*/\1/' <<<"$args"); model=$(sed -E 's/.*--model ([^ ]+).*/\1/' <<<"$args")
  seed=$(sed -E 's/.*--seed ([^ ]+).*/\1/' <<<"$args"); starts=$(sed -E 's/.*--starts ([^ ]+).*/\1/' <<<"$args")
  f=$FDABS/fits/$pair.$model.data.seed$seed.tsv
  if [ -s "$f" ] && [ $(grep -c . "$f") -ge $((starts + 1)) ]; then return 0; fi
  FD_OUT=$FDABS singularity exec --bind "${R},${scratch_dir}" $IMG python3 $CODE/fit_ext.py $args \
    > $FDABS/fits/logs/$pair.$model.seed$seed.log 2>&1
}
export -f one
mkdir -p $FDABS/fits/logs
grep . "$TASKS" | xargs -P ${NP:-64} -I{} bash -c 'one "{}"' || true
n=$(grep -c . "$TASKS"); ok=0
while read -r line; do
  pair=$(sed -E 's/.*--pair ([^ ]+).*/\1/' <<<"$line"); model=$(sed -E 's/.*--model ([^ ]+).*/\1/' <<<"$line"); seed=$(sed -E 's/.*--seed ([^ ]+).*/\1/' <<<"$line")
  f=$FDABS/fits/$pair.$model.data.seed$seed.tsv
  [ -s $f ] && [ $(grep -c . $f) -ge 6 ] && ok=$((ok+1))
done < <(grep . "$TASKS")
echo "round tasks $n, complete seed files $ok"
exit 0
