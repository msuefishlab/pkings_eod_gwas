#!/bin/bash
set -uo pipefail
source $1/env.sh
cd $O
NODE_CPUS=80 NODE_MEM=300G bash $X/submit_node.sh s3_step2 06:00:00 $X/step2_node.sh; echo "node exit $?"
cd - >/dev/null
for r in primary callable d019; do
  for f in panel.tsv ld_curves.tsv chromld.npz; do cp $O/obs/$r/$f ./obs_${r}_$f 2>/dev/null || true; done
done
cp $O/obs/primary_check.txt ./ 2>/dev/null || true
grep -iE 'error|Traceback' $O/logs/obs_rows.log | head > obs_errors.txt || true
ls -la
exit 0
