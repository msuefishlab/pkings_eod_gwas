#!/bin/bash
## submit_node.sh <job_name> <time> <script> [args] -- run <script> on one compute node and wait (called inside a
## controlling job). Account/partition names are for the MSU HPCC; edit them for your cluster. data-machine first (200 CPUs, 900G); if still PENDING after 45 min it is cancelled and the
## same script goes to general-long (64 CPUs, 240G). Job id parsed from sbatch --parsable stdout only; sacct is the
## authority on state. Returns 0 iff the node job COMPLETED.
set -uo pipefail
name=$1; tl=$2; shift 2
mkdir -p $O/logs
sub() {  # $1 = extra sbatch args, rest = script and its args
  local out extra=$1; shift
  out=$(sbatch --parsable --job-name=$name --time=$tl --output=$O/logs/${name}_%j.out $extra --wrap="bash $*" 2>>$O/logs/sbatch_stderr.txt)
  grep -oE '^[0-9]+' <<<"$out" | tail -1
}
state() { sacct -j "$1" -X -n -P -o State 2>/dev/null | head -1 | awk '{print $1}'; }
wait_job() {  # $1 jid, $2 pending limit seconds (0 = none); echo final state
  local jid=$1 lim=$2 t0=$(date +%s) st
  while true; do
    st=$(state $jid)
    case "$st" in
      COMPLETED|FAILED|CANCELLED*|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|BOOT_FAIL|DEADLINE) echo "$st"; return ;;
      PENDING|"") if [ "$lim" -gt 0 ] && [ $(( $(date +%s) - t0 )) -gt "$lim" ]; then echo "PENDING_TIMEOUT"; return; fi ;;
    esac
    sleep 60
  done
}
script=$1; shift
jid=$(sub "--account=data-machine --cpus-per-task=${NODE_CPUS:-200} --mem=${NODE_MEM:-900G}" $script "$@")
echo "[$(date)] $name data-machine job $jid"
st=$(wait_job $jid 2700)
if [ "$st" = PENDING_TIMEOUT ]; then
  scancel $jid; sleep 5
  jid=$(sub "--account=general --partition=general-long --cpus-per-task=64 --mem=240G" $script "$@")
  echo "[$(date)] $name fallback general-long job $jid"
  st=$(wait_job $jid 0)
fi
echo "[$(date)] $name job $jid final state $st"
sacct -j $jid -X -n -P -o JobID,Partition,NodeList,AllocCPUS,Elapsed,State
[ "$st" = COMPLETED ]
