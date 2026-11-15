#!/bin/bash
## 06_worker_slim.sh -- one SLiM task line from the 05* generators. The SLiM binary comes from the
## environment and the model file from $root.
## args: MODEL Q ONSET SELTPF NBPQ NTPQ MBPQ MTPQ SQ H F0 TEND SELTPQ REP SEED S_RAW
## env:  SLIM_BIN (required) -- SLiM 5.2; see README.md "Dependencies"
set -u
root="$(git rev-parse --show-toplevel)"
SLIM=${SLIM_BIN:?SLIM_BIN is not set -- point it at a SLiM 5.2 binary, see README.md}
MODEL=${SLIM_MODEL:-${root}/code/11_Drift_vs_Selection/simulation/power.slim}
line=$("$SLIM" -s "${15}" -d NBP=$5 -d NTP=$6 -d MBP=$7 -d MTP=$8 \
  -d S=$9 -d H=${10} -d F0=${11} -d SELTP=${13} -d TEND=${12} -d NSBP=67 -d NSTP=47 -d REP=${14} \
  "$MODEL" 2>/dev/null | grep '^RESULT')
if [ -n "$line" ]; then
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "${16}" "${line#RESULT	}"
fi
