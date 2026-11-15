#!/bin/bash
## 07b_shard_slim_parameterised.sh -- as 07, but the SLiM binary and model file are positional
## arguments. This is the variant used for the agreement run (run_all.sh, RUN=agree).
## args: <shard_file> <out_file> <slim_binary> <model_file>
set -u
SHARD=$1; OUT=$2; SLIM=$3; MODEL=$4
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
