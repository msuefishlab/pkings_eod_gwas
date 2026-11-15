#!/bin/bash
# orchestrate_B.sh -- Stage 2b follow-up fits, run inside ONE long controlling job, after the primary fits:
#   wide : the expanded-bound sensitivity refit (contraction <= 200x, onset 10 - 3e4 generations, same mu-union rule) of
#          each preferred model with a parameter at a bound, on BAM x TP1, same effort; nested guard also against the
#          same model at primary bounds (primary box c wide box)
#   var  : the sample-set sensitivity refits of the best alternative(s) on TP1_noAdmix (V2w_noAdm10), TP1_no017
#          (V2w_no017) and BAM_noKin x TP1 (V2w_noKin), primary bounds, same effort; spectra and the drift-test baseline
#          fits (M0-M3f) of these pairs from Stage 1 (read-only copies of $V2; waits for its FITS2.DONE). Those Stage 1
#          spectra and fits are not made by any script in this repository (see ../README.md, "Not included").
# Args: WIDE_MODELS VAR_MODELS (comma lists; "-" = none)
set -euo pipefail
WIDE_MODELS=$1; VAR_MODELS=$2
export R=${R:-${root:-$(git rev-parse --show-toplevel)}}
[ -n "${scratch_dir:-}" ] || { cd "$R" && source ./config.env && cd - >/dev/null; }
S2=${scratch_dir}/model_extension/stage2b
W=${scratch_dir}/model_extension/stage1
CODE=$R/code/11_Drift_vs_Selection/demographic_modeling/model_extension/stage2b
IMG=$R/images/dadi_tools.sif
V2=$W/build2/root/output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_var2
VPAIRS=V2w_noAdm10,V2w_no017,V2w_noKin
export CODE IMG
pids=()
if [ "$WIDE_MODELS" != "-" ]; then
  FDW=$S2/fd_wide; mkdir -p $FDW/sfs $FDW/fits
  [ -e $FDW/sfs/A1w.pkl ] || cp $S2/fd_primary/sfs/A1w.pkl $FDW/sfs/
  ls $FDW/fits/A1w.M0.data.seed1.tsv >/dev/null 2>&1 || cp $S2/fd_primary/fits/A1w.M*.data.seed*.tsv $FDW/fits/
  md5sum $FDW/sfs/A1w.pkl $S2/fd_primary/sfs/A1w.pkl > $S2/logs/wide_setup_md5.txt
  ( export RUN=wide FD=$FDW PAIRS=A1w MODELS=$WIDE_MODELS BOUNDS=wide GATE_MODELS=$WIDE_MODELS DO_LRT=0 \
           REFDIR=$S2/fd_primary/fits REFBEST=$S2/fd_primary/fits/best_fits.tsv LOGS=$S2/logs/wide
    bash $CODE/run_stage.sh > $S2/logs/orch_wide.log 2>&1 ) &
  pids+=($!)
fi
if [ "$VAR_MODELS" != "-" ]; then
  ( waited=0
    until [ -e $W/FITS2.DONE ]; do
      [ $waited -ge 2880 ] && { echo "$(date +%F_%T) no FITS2.DONE after 48 h"; exit 1; }
      sleep 600; waited=$((waited + 10))
    done
    FDV=$S2/fd_var; mkdir -p $FDV/sfs $FDV/fits
    for p in ${VPAIRS//,/ }; do
      [ -e $FDV/sfs/$p.pkl ] || cp $V2/sfs/$p.pkl $FDV/sfs/
      ls $FDV/fits/$p.M0.data.seed1.tsv >/dev/null 2>&1 || cp $V2/fits/$p.M*.data.seed*.tsv $FDV/fits/
    done
    { md5sum $FDV/sfs/*.pkl $V2/sfs/V2w_noAdm10.pkl $V2/sfs/V2w_no017.pkl $V2/sfs/V2w_noKin.pkl
      for p in ${VPAIRS//,/ }; do echo "$p baseline seed files: $(ls $FDV/fits/$p.M*.data.seed*.tsv | wc -l)"; done
      for f in $FDV/fits/V2w_*.M*.data.seed*.tsv; do cmp -s $f $V2/fits/$(basename $f) || echo "DIFF $f"; done; } > $S2/logs/var_setup_md5.txt
    export RUN=var FD=$FDV PAIRS=$VPAIRS MODELS=$VAR_MODELS BOUNDS=primary DO_LRT=0 \
           GATE_MODELS=M0,M1,M2,M3,M0f,M1f,M2f,M3f,$VAR_MODELS LOGS=$S2/logs/var
    bash $CODE/run_stage.sh > $S2/logs/orch_var.log 2>&1 ) &
  pids+=($!)
fi
rc=0
for p in "${pids[@]}"; do wait $p || rc=1; done
echo "orchestrate_B rc=$rc"
exit $rc
