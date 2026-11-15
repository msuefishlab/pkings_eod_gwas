#!/bin/bash
## orch_alt1.sh <code_dir> -- alternatives, conversion + unit check. Waits (<= 14 h) for stage2b/ADMISSIBLE.DONE, then
## alt_params.py (dadi-to-physical conversion, cross-checked against the Stage 2b physical columns); appends the
## unit-check ratio index (40) to sim_index.tsv before its first use; simulates the joint SFS of every admissible alternative and
## of baseline M2 (positive control, the drift-test agreement-run fit); dadi expected spectra (dadi_tools.sif); comparison.
set -uo pipefail
source $1/env.sh
JOBDIR=$PWD
S2=${scratch_dir}/model_extension/stage2b
IMG=$R/images/dadi_tools.sif
mkdir -p $O/alt $O/logs
w=0
until [ -e $S2/ADMISSIBLE.DONE ]; do
  [ $w -ge 840 ] && { echo "no ADMISSIBLE.DONE after 14 h"; echo "pending" > $JOBDIR/ALT_PENDING; exit 0; }
  sleep 300; w=$((w + 5))
done
echo "ADMISSIBLE.DONE seen $(date)"; cat $S2/ADMISSIBLE.DONE
cp $S2/admissible_alternatives.tsv $O/alt/admissible_alternatives.tsv
md5sum $S2/admissible_alternatives.tsv > $O/alt/admissible_alternatives.md5
$PY $X/alt_params.py $O/alt/admissible_alternatives.tsv $O/alt_params > $O/alt/alt_params.log 2>&1; echo "alt_params exit $?"
cat $O/alt/alt_params.log
cp $O/alt_params.tsv $O/alt_params.json $O/alt/alt_params.log $O/alt/admissible_alternatives.tsv $JOBDIR/ 2>/dev/null || true
MODELS=$(awk -F'\t' 'NR>1{print $1}' $O/alt_params.tsv | paste -sd, -)
if [ -z "$MODELS" ]; then echo "no admissible alternative"; touch $JOBDIR/NO_ADMISSIBLE_ALTERNATIVE; exit 0; fi
grep -q $'^ratio\t.*i_ratio 40' $O/sim_index.tsv || printf 'ratio\t\t\t\t\t\t1.0\ti_ratio 40: SFS unit check of the alternative-model conversion, physical r/mu 1.0, g 1-2 on the 100 in-sample-family blocks; appended %s before its first use; enters no criterion\t\t\n' "$(date +%F)" >> $O/sim_index.tsv
NODE_CPUS=80 NODE_MEM=300G bash $X/submit_node.sh s3_alt1 06:00:00 $X/alt1_node.sh M2,$MODELS; echo "node exit $?"
FD_OUT=$S2/fd_primary singularity exec --bind "${R},${scratch_dir}" $IMG python3 $X/alt_expected_sfs.py $R/code/11_Drift_vs_Selection/demographic_modeling/model_extension/stage2b $MODELS $O/alt/expected_alt.npz > $O/alt/expected.log 2>&1
FD_OUT=$FDA singularity exec --bind "${R},${scratch_dir}" $IMG python3 $X/alt_expected_sfs.py $R/code/11_Drift_vs_Selection/demographic_modeling/model_extension/stage2b M2 $O/alt/expected_M2.npz >> $O/alt/expected.log 2>&1
singularity exec --bind "${R},${scratch_dir}" $IMG python3 $X/alt_unit_check.py $O/alt/expected_alt.npz $O/alt/sfs $MODELS $O/alt/unit_check_alt.tsv > $O/alt/unit_check.log 2>&1
singularity exec --bind "${R},${scratch_dir}" $IMG python3 $X/alt_unit_check.py $O/alt/expected_M2.npz $O/alt/sfs M2 $O/alt/unit_check_M2.tsv >> $O/alt/unit_check.log 2>&1
cat $O/alt/expected.log $O/alt/unit_check.log
cp $O/alt/unit_check_*.tsv $O/alt/unit_check.log $O/alt/expected.log $O/sim_index.tsv $JOBDIR/ 2>/dev/null || true
grep -E 'Error|Trace' $O/alt/sfs/tasks.log | head > $JOBDIR/sfs_errors.txt || true
ls -la $JOBDIR
exit 0
