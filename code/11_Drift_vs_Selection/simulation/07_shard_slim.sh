#!/bin/bash
## 07_shard_slim.sh -- run one shard of the SLiM task file, one run per line, appending to <out_file>.
## args: <shard_file> <out_file>
## env:  SLIM_BIN (required), SLIM_MODEL (optional override of power.slim)
set -u
root="$(git rev-parse --show-toplevel)"
SHARD=$1; OUT=$2
SLIM=${SLIM_BIN:?SLIM_BIN is not set -- point it at a SLiM 5.2 binary, see README.md}
MODEL=${SLIM_MODEL:-${root}/code/11_Drift_vs_Selection/simulation/power.slim}
: > "$OUT"
while read -r M Q ONSET SELTPF NBPQ NTPQ MBPQ MTPQ SQ H F0 TEND SELTPQ REP SEED SRAW; do
  [ -z "${M:-}" ] && continue
  line=$("$SLIM" -s "$SEED" -d NBP=$NBPQ -d NTP=$NTPQ -d MBP=$MBPQ -d MTP=$MTPQ \
    -d S=$SQ -d H=$H -d F0=$F0 -d SELTP=$SELTPQ -d TEND=$TEND -d NSBP=67 -d NSTP=47 -d REP=$REP \
    "$MODEL" 2>/dev/null | grep '^RESULT')
  if [ -n "$line" ]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$M" "$Q" "$ONSET" "$SELTPF" "$SRAW" "${line#RESULT	}" >> "$OUT"
  fi
done < "$SHARD"
