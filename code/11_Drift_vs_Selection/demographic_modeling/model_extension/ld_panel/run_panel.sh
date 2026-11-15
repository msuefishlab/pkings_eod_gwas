#!/bin/bash
## run_panel.sh -- model extension, Stage 1 (LD panel): observed summary-statistic panel, one allocation.
## Usage: bash run_panel.sh <stage: main|variants>
##   main     : site table -> npz; mask-block genotype extraction (114 A1w fish); simulation/24a, 24b rerun as a check
##              (shadow root, scripts unchanged); 23a unchanged on 3 blocks (VCF identity check); panel for the samples
##              in $W/panel/<run>/samples_def.tsv; 23a/23e reproduction check.
##   variants : panel only (blocks already extracted), for $W/panel/<run>/samples_def.tsv
## Env: RUN (panel run name, default main), P (workers, default 120), PY (Python with numpy; default python3),
##      SIMOUT (outputs of code/11_Drift_vs_Selection/simulation: 23_inputs/, 24a_*.tsv, 23a_observed_ld_*.tsv,
##      23_ld_curves.tsv; default output_data/11_Drift_vs_Selection/simulation).
## Needs, from Stage 1 setup: $W/mask/{mask_blocks,chrom_split}.tsv (ext_mask_split.py) and
## $W/panel/<run>/{samples_def,pairs_def}.tsv, copied below from $PDEF (default: the shipped
## input_data/11_Drift_vs_Selection/demographic_modeling/model_extension/panel_<run>/) when absent.
set -euo pipefail
STAGE=${1:-main}
R=${root:-$(git rev-parse --show-toplevel)}
[ -n "${scratch_dir:-}" ] && [ -n "${allsites_vcf:-}" ] || { cd "$R" && source ./config.env && cd - >/dev/null; }
W=${scratch_dir}/model_extension/stage1
X=$R/code/11_Drift_vs_Selection/demographic_modeling/model_extension/ld_panel
SIM=$R/code/11_Drift_vs_Selection/simulation
FDA=$R/output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree
SIMOUT=${SIMOUT:-$R/output_data/11_Drift_vs_Selection/simulation}
VCF=${sweep_scratch}/03_QC/${allsites_vcf}
CHUNKS=$R/input_data/11_Drift_vs_Selection/demographic_modeling/chunks.tsv
PY=${PY:-python3}
RUN=${RUN:-main}; P=${P:-120}
module purge >/dev/null 2>&1; module load BCFtools/1.18-GCC-12.3.0 >/dev/null 2>&1
export LD23_DIR=$SIM OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PY X SIM FDA VCF CHUNKS W
cd $W
mkdir -p logs sitetable blocks clustering_rerun/scan ld23_check panel/$RUN/seg
## sample-set definitions: shipped for the main run (input_data/11_Drift_vs_Selection/demographic_modeling/model_extension/panel_main/); write
## panel/<run>/{samples_def,pairs_def}.tsv yourself for other runs
PDEF=${PDEF:-$R/input_data/11_Drift_vs_Selection/demographic_modeling/model_extension/panel_$RUN}
for f in samples_def.tsv pairs_def.tsv; do [ -s panel/$RUN/$f ] || [ ! -s $PDEF/$f ] || cp $PDEF/$f panel/$RUN/$f; done
which bcftools; bcftools --version | head -1
echo "start $(date)"

if [ "$STAGE" = main ]; then
  ## simulation/24b rerun (check): shadow repo root so the script's relative paths resolve read-only to the real inputs
  SH=$W/clustering_rerun/shadow
  mkdir -p $SH/output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree \
           $SH/output_data/11_Drift_vs_Selection/simulation $SH/input_data
  ln -sfn $FDA/sfs $SH/output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree/sfs
  ln -sfn $SIMOUT/23_inputs $SH/output_data/11_Drift_vs_Selection/simulation/23_inputs
  ln -sfn $R/input_data/10_Introgression $SH/input_data/10_Introgression
  ( cd $SH && $PY $SIM/24b_empirical_clustering.py $VCF $SIMOUT/24a_qualifying.tsv \
      $SIMOUT/24a_density_1kb.tsv $W/clustering_rerun/24b_clustering > $W/logs/24b_rerun.log 2>&1; \
    echo "24b exit $?" >> $W/logs/24b_rerun.log ) &
  PID24B=$!
  ## simulation/24a rerun from the on-disk site table (check; script unchanged)
  ( ls $FDA/site_table/snps.*.tsv.gz | xargs -P 8 -I{} sh -c '$PY $SIM/24a_site_scan.py {} $W/clustering_rerun/scan/$(basename {} .tsv.gz)'; \
    cat $W/clustering_rerun/scan/*.qual.tsv | sort -k1,1V -k2,2n > $W/clustering_rerun/24a_qualifying.tsv; \
    cat $W/clustering_rerun/scan/*.bins.tsv | sort -k1,1V -k2,2n > $W/clustering_rerun/24a_density_1kb.tsv; \
    md5sum $W/clustering_rerun/24a_*.tsv $SIMOUT/24a_*.tsv > $W/clustering_rerun/24a_md5.txt ) > logs/24a_rerun.log 2>&1 &
  PID24A=$!
  ## site table -> per-chunk npz
  ls $FDA/site_table/snps.*.tsv.gz | xargs -P 64 -I{} sh -c 'n=$(basename {} .tsv.gz); n=${n#snps.}; [ -s $W/sitetable/st.$n.npz ] || $PY $X/ext_sitetable.py {} $FDA/sfs/projections.tsv $FDA/sfs/deserts.tsv $W/sitetable/st.$n.npz' > logs/sitetable.log 2>&1 || echo 'WARN: sitetable tasks failed'
  echo "sitetable npz: $(ls sitetable/st.*.npz | grep -v tmp | wc -l) $(date)"
  ## union sample list = the 114 A1w fish of simulation/23_inputs
  cat $SIMOUT/23_inputs/A1w_BAM.txt $SIMOUT/23_inputs/A1w_TP1.txt > blocks/samples_union.txt
  ## genotype extraction, one task per mask block
  tail -n +2 mask/mask_blocks.tsv | awk -F'\t' '{print $1" "$2" "$5}' > blocks/tasks.txt
  grep -c . blocks/tasks.txt
  grep . blocks/tasks.txt | xargs -P $P -L 1 sh -c '[ -s $W/blocks/$0_$1.npz ] || $PY $X/ext_extract_block.py $SIM $VCF $W/blocks/samples_union.txt $CHUNKS $W/sitetable $W/blocks $0 $1 $2' > logs/extract.log 2>&1 || echo 'WARN: extraction tasks failed'
  echo "blocks npz: $(ls blocks/*.npz | grep -v tmp | wc -l) $(date)"
  ## VCF identity: 23a unchanged on three of its own blocks
  $PY $SIM/23a_observed_ld_decay.py $VCF $SIMOUT/23_inputs/A1w_BAM.txt BAM \
      $FDA/sfs/A1w.block_scan.tsv $R/input_data/10_Introgression/peaks.bed ld23_check/23a_unchanged_3blocks_BAM.tsv \
      --blocks chr1:2,chr1:5,chr25:15 > logs/23a_3blocks.log 2>&1 || echo "23a 3-block run failed"
fi

## panel
cut -f1 panel/$RUN/samples_def.tsv
ls mask/mask_blocks.tsv >/dev/null
for c in $(seq 1 25); do echo chr$c; done | xargs -P 25 -I{} sh -c '[ -s $W/panel/'$RUN'/seg/seg.{}.tsv ] && [ -s $W/panel/'$RUN'/seg/inter.{}.npz ] || $PY $X/ext_panel_chrom.py $W/blocks $W/mask/mask_blocks.tsv $W/panel/'$RUN'/samples_def.tsv $W/panel/'$RUN'/pairs_def.tsv {} $W/panel/'$RUN'/seg' > logs/panel_$RUN.log 2>&1 || echo 'WARN: panel tasks failed'
echo "seg files: $(ls panel/$RUN/seg/seg.*.tsv | wc -l) $(date)"
$PY $X/ext_aggregate.py panel/$RUN/seg mask/chrom_split.tsv panel/$RUN/samples_def.tsv panel/$RUN/pairs_def.tsv \
    panel/$RUN/panel.tsv panel/$RUN/ld_curves.tsv > logs/aggregate_$RUN.log 2>&1
if [ "$STAGE" = main ]; then
  $PY $X/ext_check_23a.py panel/$RUN/seg $SIMOUT ld23_check/check_23a.tsv > logs/check_23a.log 2>&1
  wait $PID24A || true
  wait $PID24B || true
fi
echo "done $(date)"
exit 0
