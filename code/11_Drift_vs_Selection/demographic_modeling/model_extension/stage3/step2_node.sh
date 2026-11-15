#!/bin/bash
## step2_node.sh -- Stage 3 step 2: observed rows needed by the scoring criteria. The primary all-SNV curves and bands
## are the Stage 1 panel (reused); here: 'primary' recomputed as a check (must equal Stage 1), 'callable' (LD on
## callable SNVs only, sensitivity) and 'd019' (strict-genotype mask: genotypes with DP < 8 or GQ < 20 set missing),
## each through panel_lib.segment_stats + ext_aggregate.py unchanged.
set -euo pipefail
source "$(dirname "$0")/env.sh"
NP=$(( ${SLURM_CPUS_PER_TASK:-64} - 2 ))
cd $O; mkdir -p obs logs
: > obs/tasks.txt
for row in primary callable d019; do for c in $(seq 1 25); do
  printf '%s/obs/%s/seg/inter.chr%s.npz\tobs\t%s %s %s %s %s %s chr%s %s/obs/%s/seg\n' $O $row $c $row $W/blocks $O/records $S3_MASK $O/samples_def_BAM_TP1.tsv $O/pairs_def_BAM_TP1.tsv $c $O $row >> obs/tasks.txt
done; done
grep . obs/tasks.txt | xargs -P $NP -I{} bash $X/task.sh '{}' > logs/obs_rows.log 2>&1 || echo "WARN obs tasks failed"
for row in primary callable d019; do
  $PY $X/genome_post.py $O/obs/$row $S3_SPLIT $O/samples_def_BAM_TP1.tsv $O/pairs_def_BAM_TP1.tsv >> logs/obs_rows.log 2>&1
done
## check: recomputed primary == Stage 1 panel/main (BAM, TP1, BAM|TP1 rows; LD curves)
$PY - <<'PY' > obs/primary_check.txt
import csv, os
O, W = os.environ["O"], os.environ["W"]
def rd(fn, keyf):
    return {keyf(r): r for r in csv.DictReader(open(fn), delimiter="\t")}
for name, kf, cols in (("panel.tsv", lambda r: (r["statistic"], r["sample"], r["set"]), ("value", "lo", "hi")),
                       ("ld_curves.tsv", lambda r: (r["pop"], r["set"], r["bin"]), ("r2", "lo", "hi", "n_pairs"))):
    a = rd(f"{O}/obs/primary/{name}", kf); b = rd(f"{W}/panel/main/{name}", kf)
    keys = [k for k in a if k in b]
    diff = [k for k in keys if any(a[k][c] != b[k][c] for c in cols)]
    print(name, "compared", len(keys), "of", len(a), "differing", len(diff), diff[:5])
PY
cat obs/primary_check.txt
touch STEP2.DONE
exit 0
