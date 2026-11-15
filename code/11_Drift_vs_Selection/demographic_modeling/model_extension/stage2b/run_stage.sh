#!/bin/bash
# run_stage.sh -- Stage 2b fitting driver, run inside one long controlling job. Same effort as the drift test's
# headline fits, the run_fits.sh sequence: round 1 (seeds 1-4 x 5 random starts) -> refine 2 (201-204, fold 0.5) -> nested pass 1
# (351-354) -> refine 3 (301-304, fold 0.25) -> nested pass 2 (361-364) -> collect (08 --godambe via collect_ext.py)
# -> gate_ext.py -> lrt_ext.py. Each round is one node job (sbatch, polled with sacct); maxtime 2400 s per start.
# Config (env): RUN FD PAIRS MODELS BOUNDS GATE_MODELS [REFDIR REFBEST] CODE IMG LOGS [DO_LRT]; scratch_dir and R from
# config.env (read here if not yet set). Account/partition names in submit_wait are for the MSU HPCC; edit for yours.
set -euo pipefail
: "${RUN:?}" "${FD:?}" "${PAIRS:?}" "${MODELS:?}" "${BOUNDS:?}" "${GATE_MODELS:?}" "${CODE:?}" "${IMG:?}" "${LOGS:?}"
F=$FD/fits
mkdir -p $F $LOGS
export R=${R:-${root:-$(git rev-parse --show-toplevel)}}
[ -n "${scratch_dir:-}" ] || { cd "$R" && source ./config.env && cd - >/dev/null; }
DX="singularity exec --bind ${R},${scratch_dir} $IMG python3"
REFARG=${REFDIR:+--refdir $REFDIR}

submit_wait() {  # name ncpu mem time script args...
  local name=$1 ncpu=$2 mem=$3 tlim=$4; shift 4
  local acct="--account=data-machine" fallback=0 out jid st t0
  while true; do
    out=$(sbatch --parsable $acct --job-name=$name --cpus-per-task=$ncpu --mem=$mem --time=$tlim \
          -o $LOGS/${name}_%j.log "$@" 2>>$LOGS/${name}.sbatch.err) || true
    jid=$(grep -oE '^[0-9]+' <<<"$out" | tail -1)
    [ -n "$jid" ] || { echo "$(date +%F_%T) sbatch failed for $name"; exit 1; }
    echo "$(date +%F_%T) $name submitted $jid ($acct)"
    t0=$(date +%s)
    while true; do
      sleep 60
      st=$(sacct -j "$jid" -X -n -P -o State 2>/dev/null | head -1)
      case "$st" in
        COMPLETED*) echo "$(date +%F_%T) $name $jid COMPLETED"; return 0;;
        FAILED*|CANCELLED*|TIMEOUT*|OUT_OF_ME*|NODE_FAIL*|PREEMPTED*|BOOT_FAIL*|DEADLINE*)
          echo "$(date +%F_%T) $name $jid $st"; exit 1;;
        PENDING*|"")
          if [ $fallback -eq 0 ] && [ $(( $(date +%s) - t0 )) -gt 2700 ]; then
            scancel $jid; sleep 20
            echo "$(date +%F_%T) $name $jid pending > 45 min on data-machine: cancelled, resubmitting on general-long"
            fallback=1; acct="--account=general --partition=general-long"; continue 2
          fi;;
      esac
    done
  done
}

run_round() {  # tag taskfile
  local tag=$1 tasks=$2 n
  n=$(grep -c . $tasks || true)
  echo "$(date +%F_%T) $RUN round $tag: $n tasks"
  [ "$n" -gt 0 ] || return 0
  cp $tasks $LOGS/tasks_$tag.frozen.txt
  local np=$(( n < 48 ? n : 48 ))
  submit_wait s2b_${RUN}_$tag $np $((np * 4))G 10:00:00 --wrap \
    "NP=$np bash $CODE/fit_round_ext.sh $LOGS/tasks_$tag.frozen.txt $FD $CODE $IMG"
}

gen() {  # mode n out record
  $DX $CODE/tasks_ext.py $1 --n $2 --fits $F --sfs $FD/sfs --pairs $PAIRS --models $MODELS --bounds $BOUNDS \
     $REFARG --out $3 --record $4
}

if [ ! -e $F/ROUND1.DONE ]; then
  gen round1 0 $F/tasks_fit.txt $F/round1_thetaref.txt
  run_round r1 $F/tasks_fit.txt; touch $F/ROUND1.DONE
fi
refine() {
  [ -e $F/ROUND$1.DONE ] && return 0
  gen refine $1 $F/tasks_refine$1.txt $F/refine$1_starts.txt
  run_round r$1 $F/tasks_refine$1.txt; touch $F/ROUND$1.DONE
}
nested() {
  [ -e $F/NESTED$1.DONE ] && return 0
  gen nested $1 $F/tasks_nested$1.txt $F/nested$1_starts.txt
  echo "nested pass $1 violations:"; cat $F/nested$1_starts.txt
  run_round n$1 $F/tasks_nested$1.txt; touch $F/NESTED$1.DONE
}
refine 2
nested 1
refine 3
nested 2

if [ ! -e $F/COLLECT.DONE ]; then
  cat > $LOGS/collect_$RUN.sh <<CS
#!/bin/bash
set -euo pipefail
cd $(dirname $FD); export OMP_NUM_THREADS=1
FD_OUT=$FD $DX $CODE/collect_ext.py --godambe
FD_OUT=$FD $DX $CODE/gate_ext.py --pairs $PAIRS --models $GATE_MODELS --bounds $BOUNDS ${REFBEST:+--refdir $REFBEST} --out $FD/gate_ext.tsv
if [ "${DO_LRT:-1}" = "1" ]; then FD_OUT=$FD $DX $CODE/lrt_ext.py --pairs $PAIRS --out $FD/lrt_ext.tsv; fi
CS
  submit_wait s2b_${RUN}_collect 24 192G 16:00:00 $LOGS/collect_$RUN.sh
  touch $F/COLLECT.DONE
fi
for t in fit refine2 refine3 nested1 nested2; do [ -e $F/tasks_$t.txt ] && echo "$t $(grep -c . $F/tasks_$t.txt || true)"; done > $FD/task_counts.txt
echo "seed files $(ls $F/*.data.seed*.tsv | wc -l)" >> $FD/task_counts.txt
touch $FD/STAGE.DONE
echo "$(date +%F_%T) $RUN done"
