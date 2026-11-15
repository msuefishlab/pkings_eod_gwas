#!/bin/bash
# orchestrate.sh — smc++ Ne(t) + split-time pipeline on callability-corrected masks
# (masks_v2), 8 populations. One controlling job; SLURM arrays are submitted inside.
#
#   Stage 1   vcf2smc, single population   (8 pops x 25 chr x 5 DIs = 1000 tasks)
#   Stage 5a  vcf2smc, joint (split input) (6 pairs x 25 chr x 2 DIs = 300 tasks)  [concurrent]
#   Stage 3   smc++ estimate               (8 pops x 2 mu = 16)  -> memory gate
#   Stage 4   bootstrap estimates          (10 x 16 = 160)                         [concurrent]
#   Stage 5b  smc++ split                  (6 pairs x 2 mu x 2 replicate runs = 24)
#   Collate   results/*.tsv + per-pop Ne_tables + resource_usage.tsv -> $OUTD
#
# Usage (from the repository root, after `source config.env`):
#   sbatch code/12_Historical_Ne/scripts/run_smcpp.sb
# Idempotent: every stage rebuilds its to-do list from what is already on disk and
# adopts a still-running array of the same name, so resubmitting resumes.
# DRYRUN=1 bash orchestrate.sh  -> pre-flight + manifests only, nothing submitted.
#
# Env overrides: REPO, SMCPP_WORK (smc++ working tree, default ${scratch_dir}/smcpp),
# V2 (run tree), OUTD (export dir), PY (python with numpy),
# CB (dir with bedtools/bgzip/tabix; default: the directory of bedtools on PATH).
set -uo pipefail
SCR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)       # code/12_Historical_Ne/scripts
REPO=${REPO:-$(cd "$SCR/../../.." && pwd)}
[ -n "${scratch_dir:-}" ] || { pushd "$REPO" >/dev/null; source config.env; popd >/dev/null; }
export SMCPP_WORK=${SMCPP_WORK:-${scratch_dir}/smcpp}   # read by the w_*.sh array workers
W=$SMCPP_WORK
INP=$REPO/input_data/12_Historical_Ne    # sample lists and the vcf2smc task manifest
SIF=$REPO/images/smcpp.sif
V2=${V2:-$W/v2}
DRY=${DRYRUN:-0}
LOGS=$V2/logs; MAN=$V2/manifests
PY=${PY:-python3}
CB=${CB:-$(dirname "$(command -v bedtools)")}
OUTD=${OUTD:-$V2/export}
POPS="BAM BIK BAVA BIR APA_2009 BEN DOV DOG"
PAIRS="BAM:APA_2009 BIK:DOV BAVA:DOV BIR:DOV APA_2009:BEN DOV:DOG"
mkdir -p "$LOGS" "$MAN" "$V2/smc_input" "$V2/estimates" "$V2/results"
for v in $(env | grep -o '^SBATCH_[A-Z_]*' || true); do unset "$v"; done
unset SLURM_JOB_PARTITION || true
log() { echo "[$(date "+%F %T")] $*" | tee -a "$V2/STATUS.log" >&2; }
nlines() { grep -c . "$1" 2>/dev/null || true; }

# ------------------------------------------------------------ array runner --
LAST_JID=""
wait_array() {  # jid name n
  local jid=$1 name=$2 n=$3 st="" act i=0
  sleep 60
  while :; do
    st=$(sacct -X -n -P -j "$jid" -o State 2>/dev/null || true)
    if [ -n "$st" ]; then
      act=$(grep -cE 'PENDING|RUNNING|REQUEUED|CONFIGURING|COMPLETING|RESIZING|SUSPENDED' <<<"$st" || true)
      [ "$act" -eq 0 ] && break
      [ $((i % 10)) -eq 0 ] && log "$name: $(grep -c '^COMPLETED' <<<"$st" || true)/$n completed, $act active"
    fi
    i=$((i + 1)); sleep 90
  done
  log "$name: array $jid finished — $(sed 's/ by .*//' <<<"$st" | sort | uniq -c | tr -s ' ' | tr '\n' ';')"
}
find_active() {  # print base job id of an active array named v2_<name>, if any
  sacct -u "$USER" -X -n -P -S now-4days --name="v2_$1" -o JobID,State 2>/dev/null \
    | grep -E 'PENDING|RUNNING|REQUEUED|CONFIGURING|COMPLETING|SUSPENDED' \
    | head -1 | cut -d'|' -f1 | grep -oE '^[0-9]+' || true
}
run_array() {   # name manifest worker conc partition time mem cpus
  local name=$1 man=$2 wk=$3 conc=$4 part=$5 tl=$6 mem=$7 cpu=$8 n out ex frozen tries=0
  LAST_JID=""
  if [ "$DRY" != 1 ]; then
    ex=$(find_active "$name")
    if [ -n "$ex" ]; then   # resume: an array for this stage is already in flight
      LAST_JID=$ex; grep -q "^$name	$ex	" "$V2/jobids.tsv" 2>/dev/null || echo -e "$name\t$ex\tadopted" >> "$V2/jobids.tsv"
      log "$name: adopting active array $ex (not resubmitting)"
      wait_array "$ex" "$name" "?"; return 0
    fi
  fi
  n=$(nlines "$man")
  if [ "$n" -eq 0 ]; then log "$name: nothing to do"; return 0; fi
  if [ "$DRY" = 1 ]; then log "DRY $name: $n tasks %$conc $part $tl $mem ${cpu}c"; return 0; fi
  mkdir -p "$LOGS/slurm_$name"
  # freeze the manifest: running tasks read lines by index, so the file an
  # array was submitted with must never be rewritten afterwards
  frozen="$MAN/frozen_${name}_$(date +%Y%m%d%H%M%S).tsv"; cp "$man" "$frozen"
  while :; do   # retry covers MaxSubmitJobs=1000 while another array is queued
    out=$(sbatch --parsable --array=1-${n}%${conc} --account=general --partition="$part" \
          --time="$tl" --mem="$mem" --cpus-per-task="$cpu" -J "v2_$name" \
          -o "$LOGS/slurm_$name/%A_%a.out" "$wk" "$frozen" "$V2" "$REPO" 2>"$LOGS/sbatch_$name.err")
    LAST_JID=$(grep -oE '^[0-9]+' <<<"$out" | tail -1)
    [ -n "$LAST_JID" ] && break
    tries=$((tries + 1))
    { [ $tries -le 2 ] || [ $((tries % 15)) -eq 0 ]; } && log "$name: sbatch refused (try $tries): $(grep -v -e Replacing -e 'more information' -e 'docs.icer' "$LOGS/sbatch_$name.err" | tail -c 200)"
    [ $tries -ge 720 ] && { log "$name: GAVE UP submitting"; return 1; }
    sleep 120
  done
  echo -e "$name\t$LAST_JID\t$n" >> "$V2/jobids.tsv"
  log "$name: array $LAST_JID submitted ($n tasks, %$conc, $part, $tl, $mem, ${cpu}c, manifest $(basename "$frozen"))"
  wait_array "$LAST_JID" "$name" "$n"
}

# ----------------------------------------------------------- to-do builders --
todo_main() {  # full LPT manifest minus outputs on disk
  tr ' ' '\t' < "$INP/vcf2smc_tasks_sub20_lpt.txt" | grep . | while IFS=$'\t' read -r P C D I; do
    [ -s "$V2/smc_input/$P/${C}_di${I}.smc.gz" ] || printf '%s\t%s\t%s\t%s\n' "$P" "$C" "$D" "$I"; done > "$MAN/todo_main.tsv"
  nlines "$MAN/todo_main.tsv"
}
todo_joint() {
  local chrs; chrs=$(awk '{print $2}' "$INP/vcf2smc_tasks_sub20_lpt.txt" | awk '!s[$0]++')
  : > "$MAN/todo_joint.tsv"
  for c in $chrs; do for pd in $PAIRS; do
    p1=${pd%%:*}; p2=${pd##*:}; pr="$p1-$p2"
    d1=$(head -1 "$INP/samples/${p1}_sub20.txt"); d2=$(head -1 "$INP/samples/${p2}_sub20.txt")
    [ -s "$V2/smc_input_split/$pr/${c}_pop1di.smc.gz" ] || printf '%s\t%s\t%s\t%s\tpop1di\t%s\n' "$pr" "$p1" "$p2" "$c" "$d1" >> "$MAN/todo_joint.tsv"
    [ -s "$V2/smc_input_split/$pr/${c}_pop2di.smc.gz" ] || printf '%s\t%s\t%s\t%s\tpop2di\t%s\n' "$pr" "$p1" "$p2" "$c" "$d2" >> "$MAN/todo_joint.tsv"
  done; done
  nlines "$MAN/todo_joint.tsv"
}
pop_complete() { [ "$(ls "$V2/smc_input/$1"/chr*_di*.smc.gz 2>/dev/null | wc -l)" -eq 125 ]; }
todo_point() {
  : > "$MAN/todo_point.tsv"
  for p in $POPS; do pop_complete "$p" || { log "estimates: $p marginal inputs incomplete — skipped"; continue; }
    for a in "1e-8 1e8" "4e-9 4e9"; do set -- $a
      [ -s "$V2/estimates/${p}_mu$2/model.final.json" ] || printf 'point\t0\t%s\t%s\t%s\n' "$p" "$1" "$2" >> "$MAN/todo_point.tsv"
    done; done
  nlines "$MAN/todo_point.tsv"
}
todo_boot() {
  grep . "$MAN/boot_all.tsv" 2>/dev/null | while IFS=$'\t' read -r k K p mu lab; do
    [ -s "$V2/bootstraps/estimates/boot_$K/${p}_mu$lab/model.final.json" ] || printf '%s\t%s\t%s\t%s\t%s\n' "$k" "$K" "$p" "$mu" "$lab"; done > "$MAN/todo_boot.tsv"
  nlines "$MAN/todo_boot.tsv"
}
todo_split() {
  : > "$MAN/todo_split.tsv"
  for pd in $PAIRS; do p1=${pd%%:*}; p2=${pd##*:}; pr="$p1-$p2"
    for mu in 1e8 4e9; do for mode in unfolded folded; do
      sub=split; [ $mode = folded ] && sub=split_folded
      [ -s "$V2/$sub/${pr}_mu$mu/model.final.json" ] && continue
      [ -s "$V2/estimates/${p1}_mu$mu/model.final.json" ] && [ -s "$V2/estimates/${p2}_mu$mu/model.final.json" ] || continue
      [ "$(ls "$V2/smc_input_split/$pr"/*.smc.gz 2>/dev/null | wc -l)" -eq 50 ] || continue
      printf '%s\t%s\t%s\tmu%s\t%s\n' "$pr" "$p1" "$p2" "$mu" "$mode" >> "$MAN/todo_split.tsv"
    done; done; done
  nlines "$MAN/todo_split.tsv"
}

# ================================================================ STAGE 0 ====
if [ -f "$V2/jobids.tsv" ]; then grep -E '^[a-z_]+	[0-9]+	' "$V2/jobids.tsv" > "$V2/jobids.tsv.clean" || true; mv "$V2/jobids.tsv.clean" "$V2/jobids.tsv"; fi
log "==== v2 rerun start (job ${SLURM_JOB_ID:-none}, DRY=$DRY) V2=$V2"
fail=0
[ -s "$SIF" ] || { log "PREFLIGHT missing $SIF"; fail=1; }
for p in $POPS; do
  [ -s "$W/audit/masks_v2/$p.mask.bed.gz" ] && [ -s "$W/audit/masks_v2/$p.mask.bed.gz.tbi" ] || { log "PREFLIGHT missing mask $p"; fail=1; }
  [ -s "$INP/samples/${p}_sub20.txt" ] || { log "PREFLIGHT missing samples $p"; fail=1; }
done
nv=$(ls "$W"/data/vcf_by_chr/chr*.vcf.gz 2>/dev/null | wc -l); nt=$(ls "$W"/data/vcf_by_chr/chr*.vcf.gz.tbi 2>/dev/null | wc -l)
[ "$nv" -eq 25 ] && [ "$nt" -eq 25 ] || { log "PREFLIGHT per-chr VCFs $nv/25, tbi $nt/25"; fail=1; }
[ "$(nlines "$INP/vcf2smc_tasks_sub20_lpt.txt")" -eq 1000 ] || { log "PREFLIGHT manifest not 1000 lines"; fail=1; }
# Optional: reuse BAM's 125 masks_v2 .smc.gz files from an earlier run with the identical
# command (audit/rerun_v2/smc_input/BAM) if present.
nb=$(ls "$W"/audit/rerun_v2/smc_input/BAM/*.smc.gz 2>/dev/null | wc -l)
log "preflight: masks 8, VCFs $nv/25, BAM pilot files $nb/125, scratch free $(df -h "$W" | awk 'NR==2{print $4}')"
[ $fail -eq 0 ] || { log "PREFLIGHT FAILED — abort"; exit 2; }
mkdir -p "$V2/smc_input/BAM"
if [ "$DRY" = 1 ] || [ "$nb" -ne 125 ]; then log "BAM pilot staging skipped (DRY=$DRY, pilot files $nb/125)"
else
  for f in "$W"/audit/rerun_v2/smc_input/BAM/*.smc.gz; do
    d="$V2/smc_input/BAM/$(basename "$f")"
    [ -s "$d" ] || { cp -p "$f" "$d.tmp.smc.gz" && gzip -t "$d.tmp.smc.gz" && mv "$d.tmp.smc.gz" "$d"; }
  done
  log "staged BAM: $(ls "$V2"/smc_input/BAM/*.smc.gz | wc -l)/125"
fi
log "to-do: main vcf2smc $(todo_main)/1000, joint vcf2smc $(todo_joint)/300"

# ============================================ STAGE 1 & 5a (concurrent) =====
stage_main() {
  run_array main "$MAN/todo_main.tsv" "$SCR/w_vcf2smc.sh" 500 general-short 03:50:00 8G 1
  [ "$DRY" = 1 ] && return 0
  n=$(todo_main); [ "$n" -gt 0 ] && { log "main: $n missing → retry on general-long"; \
    run_array main_retry "$MAN/todo_main.tsv" "$SCR/w_vcf2smc.sh" 500 general-long 12:00:00 16G 1; }
  log "STAGE 1 done: main vcf2smc missing=$(todo_main)/1000"
}
build_masks() {
  mkdir -p "$V2/masks_split"
  for pd in $PAIRS; do p1=${pd%%:*}; p2=${pd##*:}; m="$V2/masks_split/${p1}_${p2}_union.mask.bed.gz"
    [ -s "$m" ] && [ -s "$m.tbi" ] && continue
    [ "$DRY" = 1 ] && { log "DRY: would build $m"; continue; }
    { zcat "$W/audit/masks_v2/$p1.mask.bed.gz"; zcat "$W/audit/masks_v2/$p2.mask.bed.gz"; } \
      | LC_ALL=C sort -k1,1 -k2,2n -S 6G --parallel=4 -T "${TMPDIR:-$V2}" \
      | "$CB/bedtools" merge -i - | "$CB/bgzip" > "$m.tmp" && mv "$m.tmp" "$m" && "$CB/tabix" -f -p bed "$m"
    log "mask ${p1}_${p2}: $(zcat "$m" | wc -l) regions"
  done
}
stage_joint() {
  build_masks
  todo_joint >/dev/null
  run_array joint "$MAN/todo_joint.tsv" "$SCR/w_joint.sh" 300 general-long 12:00:00 12G 1
  [ "$DRY" = 1 ] && return 0
  n=$(todo_joint); [ "$n" -gt 0 ] && { log "joint: $n missing → retry"; \
    run_array joint_retry "$MAN/todo_joint.tsv" "$SCR/w_joint.sh" 300 general-long 24:00:00 24G 1; }
  log "STAGE 5a done: joint vcf2smc missing=$(todo_joint)/300"
}
todo_main >/dev/null
stage_main & PM=$!
sleep 5
stage_joint & PJ=$!
wait $PM

# ========================================= STAGE 2+3: estimates (mem gate) ==
n=$(todo_point); log "STAGE 3: $n point estimates to run"
run_array estimate "$MAN/todo_point.tsv" "$SCR/w_estimate.sh" 16 general-short 04:00:00 64G 4
EJ=$LAST_JID
n=$(todo_point); [ "$n" -gt 0 ] && [ "$DRY" != 1 ] && { log "estimate: $n missing → retry with 192G"; \
  run_array estimate_retry "$MAN/todo_point.tsv" "$SCR/w_estimate.sh" 16 general-long 24:00:00 192G 4; }
# Memory gate: size bootstraps from the measured peak RSS of the 16 point estimates
# (default 39 GB, the peak RSS of a BAM masks_v2 estimate).
maxg=39
if [ -n "$EJ" ]; then
  m=$(sacct -j "$EJ" -n -P --units=G -o MaxRSS 2>/dev/null | sed 's/G$//' | awk '$1+0>x{x=$1+0}END{printf "%d", x+0.999}')
  [ -n "$m" ] && [ "$m" -gt 0 ] && maxg=$m
fi
BOOTMEM=$(( maxg * 3 / 2 )); [ $BOOTMEM -lt 64 ] && BOOTMEM=64; [ $BOOTMEM -gt 480 ] && BOOTMEM=480
log "MEMORY GATE: max estimate RSS ≈ ${maxg} GB → bootstrap request ${BOOTMEM}G per task"
log "STAGE 3 done: point estimates missing=$(todo_point)"

# ======================================== STAGE 4 & 5b (concurrent) ========
stage_boot() {
  [ "$DRY" = 1 ] && { log "DRY: bootstraps would be built after stage 1"; return 0; }
  "$PY" "$SCR/make_bootstraps.py" "$V2" "$MAN/boot_all.tsv" 2>&1 | tee -a "$V2/STATUS.log"
  todo_boot >/dev/null
  run_array boot "$MAN/todo_boot.tsv" "$SCR/w_estimate.sh" 160 general-short 04:00:00 ${BOOTMEM}G 4
  n=$(todo_boot); [ "$n" -gt 0 ] && { log "boot: $n missing → retry"; \
    run_array boot_retry "$MAN/todo_boot.tsv" "$SCR/w_estimate.sh" 160 general-long 24:00:00 $((BOOTMEM*2))G 4; }
  log "STAGE 4 done: bootstraps missing=$(todo_boot)/160"
}
stage_boot & PB=$!
wait $PJ
n=$(todo_split); log "STAGE 5b: $n split runs (unfolded + folded)"
run_array split "$MAN/todo_split.tsv" "$SCR/w_split.sh" 24 general-long 24:00:00 96G 8
n=$(todo_split); [ "$n" -gt 0 ] && [ "$DRY" != 1 ] && { log "split: $n missing → retry"; \
  run_array split_retry "$MAN/todo_split.tsv" "$SCR/w_split.sh" 24 general-long 72:00:00 256G 8; }
log "STAGE 5b done: split runs missing=$(todo_split)"
wait $PB

# ================================================================ COLLATE ===
[ "$DRY" = 1 ] && { log "DRY run complete"; exit 0; }
"$PY" "$SCR/collate.py" "$V2" "$V2/results" 2>&1 | tee -a "$V2/STATUS.log"
{ echo -e "stage\tJobID\tJobName\tPartition\tAllocCPUS\tReqMem\tMaxRSS\tElapsed\tState"
  grep . "$V2/jobids.tsv" 2>/dev/null | while IFS=$'\t' read -r nm j n; do
    sacct -j "$j" -n -P --units=G -o JobID,JobName,Partition,AllocCPUS,ReqMem,MaxRSS,Elapsed,State 2>/dev/null | sed "s/^/$nm\t/" | tr '|' '\t'
  done; } > "$V2/results/resource_usage.tsv"
mkdir -p "$OUTD"
mkdir -p "$OUTD/Ne_tables"; cp "$V2"/results/Ne_tables/*.tsv "$OUTD/Ne_tables/" 2>/dev/null || true
cp "$V2"/results/*.tsv "$V2"/results/collation_status.txt "$V2/STATUS.log" "$V2/jobids.tsv" "$OUTD"/ 2>/dev/null || true
for d in "$V2"/estimates/*/; do b=$(basename "$d"); cp "$d/model.final.json" "$OUTD/est__$b.json" 2>/dev/null || true; done
for s in split split_folded; do for d in "$V2/$s"/*/; do [ -d "$d" ] || continue; b=$(basename "$d"); cp "$d/model.final.json" "$OUTD/${s}__$b.json" 2>/dev/null || true; done; done
for d in "$V2"/bootstraps/estimates/boot_*/*/; do k=$(basename "$(dirname "$d")"); b=$(basename "$d"); cp "$d/model.final.json" "$OUTD/boot__${k}__$b.json" 2>/dev/null || true; done
log "exported $(ls "$OUTD" | wc -l) files to $OUTD"
if grep -q 'ALL PRESENT' "$V2/results/collation_status.txt"; then log "==== COMPLETE: all outputs present"; exit 0
else log "==== INCOMPLETE: $(grep -c . "$V2/results/collation_status.txt") items missing (see collation_status.txt); resubmit to resume"; exit 3; fi
