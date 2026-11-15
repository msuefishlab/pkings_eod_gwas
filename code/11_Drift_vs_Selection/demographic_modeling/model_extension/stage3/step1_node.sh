#!/bin/bash
## step1_node.sh -- Stage 3 step 1 (pipeline), on one node: per-record site tables of the 860 M-sum blocks
## (records_full.py), their checks, the forward heterozygote-dropout table, and the coord-pilot b'' reproduction
## under M2 (a check against the pilot run's own simulations).
set -euo pipefail
source "$(dirname "$0")/env.sh"
module purge >/dev/null 2>&1; module load BCFtools/1.18-GCC-12.3.0 >/dev/null 2>&1
NP=$(( ${SLURM_CPUS_PER_TASK:-64} - 4 ))
cd $O; mkdir -p records logs pilot_repro
{ echo "host $(hostname)"; lscpu | grep -E 'Model name|^CPU\(s\)'; echo "cpus ${SLURM_CPUS_PER_TASK:-} start $(date)"
  bcftools --version | head -1
  $PY -c 'import msprime, numpy, tskit; print("msprime", msprime.__version__, "numpy", numpy.__version__, "tskit", tskit.__version__)'
  md5sum $PANEL_DIR/panel_lib.py $PANEL_DIR/ext_aggregate.py $LD23_DIR/ld23_lib.py $LD23_DIR/23e_fit_r_over_mu.py $S3_PARAMS $X/*.py $X/*.sh; } > logs/provenance_step1.txt 2>&1
## 1a per-record site tables
tail -n +2 $S3_MASK | awk -F'\t' -v O=$O -v A="$VCF $INP/groups.tsv $INP/group_sizes.tsv $W/blocks/samples_union.txt $FDA/sfs/deserts.tsv $W/blocks $PILOT/a2/records" \
  '{printf "%s/records/%s_%s.rec.npz\trec\t%s %s %s %s %s/records\n", O, $1, $2, A, $1, $2, $6, O}' > records/tasks.txt
echo "record tasks $(grep -c . records/tasks.txt) $(date)"
grep . records/tasks.txt | xargs -P $NP -I{} bash $X/task.sh '{}' > logs/records.log 2>&1 || echo "WARN record tasks failed"
$PY - <<'PY'
import csv, glob, os
O = os.environ["O"]
rows = []
for f in sorted(glob.glob(f"{O}/records/*.check.tsv")):
    rows.append(dict(zip(*[l.rstrip("\n").split("\t") for l in open(f)])))
keys = list(rows[0]); [keys.append(k) for r in rows for k in r if k not in keys]
with open(f"{O}/records_checks.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=keys, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
ok = [r for r in rows if r["L_match"] == "1" and r["npz_snps"] == r["npz_snps_found"] == r["callable_flag_match"]
      and r["geno_compared"] == r["geno_miss_match"] and r["npz_snps_found"] == r["npz_snps_have_gq"]
      and r["pilot_compared"] in ("0", "1") and (r["pilot_compared"] == "0" or
      (r["pilot_rpos_identical"] == r["pilot_miss_identical"] == r["pilot_callable_identical"] == "1"))]
pil = sum(r["pilot_compared"] == "1" for r in ok)
L = sum(int(r["n_callable_unique"]) for r in rows)
print(f"record blocks {len(rows)}, passing {len(ok)}, pilot blocks identical {pil}, callable bp {L}")
open(f"{O}/records_summary.txt", "w").write(f"blocks {len(rows)} passing {len(ok)} pilot_identical {pil} callable_bp {L}\n")
if len(ok) != 860 or pil != 40 or L != 606107646:
    open(f"{O}/RECORDS_MISMATCH", "w").write("see records_checks.tsv\n")
PY
if [ -e RECORDS_MISMATCH ]; then echo "STOP: record checks failed $(date)"; exit 0; fi
## 1b forward heterozygote-dropout table
head -2 $W/panel/main/samples_def.tsv > samples_def_BAM_TP1.tsv
printf 'BAM\tTP1\n' > pairs_def_BAM_TP1.tsv
[ -s dropout_d.tsv ] || $PY $X/dropout_rates.py $W/blocks $O/records $S3_MASK $W/blocks/samples_union.txt samples_def_BAM_TP1.tsv dropout_d.tsv 64 > logs/dropout.log 2>&1
echo "dropout table: $(grep -c . dropout_d.tsv) lines $(date)"
## 1c records cache + pilot reproduction (M2, pilot seeds, 40 blocks x 5 r/mu x 100 reps)
free=$(df -BG --output=avail ${TMPDIR:-/tmp} | tail -1 | tr -dc 0-9)
if [ "${free:-0}" -gt 150 ]; then CACHE=${TMPDIR}/s3cache; else CACHE=$O/cache; fi
$PY $X/prep_cache.py $O/records $S3_MASK $CACHE 64 > logs/cache_step1.log 2>&1
export S3_CACHE=$CACHE
PB=$R/code/11_Drift_vs_Selection/demographic_modeling/model_extension/coord_pilot/pilot_blocks.tsv
: > pilot_repro/tasks.txt
for ir in 4 3 2 1 0; do for ib in $(seq 0 39); do
  printf '%s/pilot_repro/pilot_%s_%s.npz\tsim\tpilot %s %s %s 0.0237707 %s/pilot_repro/pilot_%s_%s.npz\n' $O $ib $ir $ib $ir $PB $O $ib $ir >> pilot_repro/tasks.txt
done; done
echo "pilot tasks $(grep -c . pilot_repro/tasks.txt) $(date)"
grep . pilot_repro/tasks.txt | xargs -P $NP -I{} bash $X/task.sh '{}' > logs/pilot_repro.log 2>&1 || echo "WARN pilot tasks failed"
$PY $X/pilot_compare.py $O/pilot_repro $PILOT/a2/sim $W/ld_track/observed_ld_curves.tsv $O/pilot_repro.tsv > logs/pilot_compare.log 2>&1
cat logs/pilot_compare.log
## 1d whole-genome path test (M2, physical r/mu 1.4, i_ratio 99, g 1): 25 chromosome tasks + genome_post (format/timing)
TG=$O/test_genome/M2/r99/g1
for c in $(seq 1 25); do printf '%s/seg/inter.chr%s.npz\tsim\tgenome M2 99 1.4 1 chr%s %s/seg\n' $TG $c $c $TG; done > logs/test_genome_tasks.txt
grep . logs/test_genome_tasks.txt | xargs -P 25 -I{} bash $X/task.sh '{}' > logs/test_genome.log 2>&1 || echo "WARN test genome failed"
$PY $X/genome_post.py $TG $S3_SPLIT samples_def_BAM_TP1.tsv pairs_def_BAM_TP1.tsv >> logs/test_genome.log 2>&1 || echo "WARN test post failed"
grep -E 'WALL|cpu=' logs/test_genome.log | sort -k3 -n | tail -3
touch STEP1.DONE
echo "done $(date)"
exit 0
