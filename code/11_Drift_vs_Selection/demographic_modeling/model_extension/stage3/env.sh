## env.sh -- Stage 3 paths (sourced by every Stage 3 shell script). Reads config.env from the repository root if it
## has not been sourced yet. PY must be a Python with msprime 1.4.4, tskit and numpy (default: python3 on PATH).
export R=${R:-${root:-$(git rev-parse --show-toplevel)}}
[ -n "${scratch_dir:-}" ] && [ -n "${allsites_vcf:-}" ] || { cd "$R" && source ./config.env && cd - >/dev/null; }
export W=${scratch_dir}/model_extension/stage1
export O=${scratch_dir}/model_extension/stage3
export X=$R/code/11_Drift_vs_Selection/demographic_modeling/model_extension/stage3
export PILOT=${scratch_dir}/model_extension/coord_pilot
export PANEL_DIR=$R/code/11_Drift_vs_Selection/demographic_modeling/model_extension/ld_panel
export LD23_DIR=$R/code/11_Drift_vs_Selection/simulation
export SIM=$LD23_DIR
export FDA=$R/output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree
export INP=$R/input_data/11_Drift_vs_Selection/demographic_modeling
export VCF=${sweep_scratch}/03_QC/${allsites_vcf}
export PY=${PY:-python3}
export S3_MASK=$W/mask/mask_blocks.tsv S3_SPLIT=$W/mask/chrom_split.tsv S3_INDEX=$O/sim_index.tsv
export S3_PARAMS=$SIM/params.as_run.agree.json S3_DROPOUT=$O/dropout_d.tsv S3_ALT=$O/alt_params.tsv
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
