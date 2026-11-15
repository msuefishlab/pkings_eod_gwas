# model_extension: contraction histories and linked simulations for the Bambomo drift test

Produces: the ESM "Drift Test at GWAS Peaks" model-based regional test (exploratory in the ESM): two-population
histories that allow population contraction, fitted to the Bambomo x TP1 (Apassa + Bengue) joint spectrum, and
simulations of linked variation under the fitted histories, scored against observed linkage disequilibrium (LD) at
short distances and out to 100 kb. These simulations did not reproduce the observed LD decay in Bambomo, so no neutral
probability is assigned to the Peak 1 region.

- `stage2b/`: six alternative models, three contraction families each fitted to the unfolded and the folded spectrum
  (`CB`: BP size change at T_c; `CBT`: BP and TP change; `CBTM`: `CBT` plus asymmetric migration that stops at T_m),
  fitted next to the drift test's eight models (M0-M3, M0f-M3f). Fitting rounds, convergence, the adequacy gate,
  Godambe-adjusted LRTs (information only) and the tables.
- `ld_panel/`: the observed summary-statistic panel (pi, dxy, Fst, Delta ladder, binned r^2 to 100 kb,
  interchromosomal r^2) on the callable 1-Mb blocks, with chromosome-bootstrap bands and a calibration / held-out
  chromosome split.
- `stage3/`: msprime simulations of every mask block under each history, r/mu fits against the calibration
  chromosomes (criterion P1: held-out LD curve inside the observed band in >= 80% of 2-kb bins; P3: held-out mean r^2
  over 20-100 kb inside the simulated central 95%), the false-rejection-rate (FRR) calibration of these criteria on
  pseudo-observed genomes, and the scoring of the real data.
- `coord_pilot/pilot_blocks.tsv`: the 40 blocks of the coordinate pilot, read by Stage 3 step 1 for its reproduction
  check.

Software: `images/dadi_tools.sif` (dadi; Stage 2b and the SFS unit check), a Python with msprime 1.4.4, tskit, numpy,
pandas and matplotlib (`PY`, default `python3`), bcftools 1.18. Shell drivers submit SLURM jobs; edit the account and
partition names in `stage2b/run_stage.sh` and `stage3/submit_node.sh` for your cluster.

## Inputs
- `output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree/` (`$FDA`): the drift test run with
  `FD_OUT` set to this directory (`../fixed_differences/`): `site_table/`, `sfs/` (`A1w.pkl`, `A1w.block_scan.tsv`,
  `projections.tsv`, `deserts.tsv`), `fits/` (seed files, `best_fits.tsv`).
- `input_data/11_Drift_vs_Selection/demographic_modeling/{chunks,groups,group_sizes}.tsv`: from
  `../fixed_differences/00_prepare_inputs.py`.
- `code/11_Drift_vs_Selection/simulation/`: `23a_observed_ld_decay.py`, `23e_fit_r_over_mu.py`, `ld23_lib.py`,
  `24a_site_scan.py`, `24b_empirical_clustering.py`, `params.as_run.agree.json`; and that module's outputs (`SIMOUT`,
  default `output_data/11_Drift_vs_Selection/simulation/`): `23_inputs/A1w_{BAM,TP1}.txt`, `24a_*.tsv`,
  `23a_observed_ld_{BAM,TP1}.tsv`, `23_ld_curves.tsv`.
- `${sweep_scratch}/03_QC/${allsites_vcf}`: the all-sites VCF, staged by `sweep_detection/00_fetch_inputs.sb`.
- `input_data/10_Introgression/peaks.bed` and the reference `.fai` (`input_data/00_Reference_Genome/`).

Working trees: `${scratch_dir}/model_extension/{stage1,stage2b,stage3,coord_pilot}` (`W`, `S2`, `O`, `PILOT`).

## Steps
Run from the repository root after `source config.env`; `R=$root`, `W=${scratch_dir}/model_extension/stage1`,
`S2=${scratch_dir}/model_extension/stage2b`, `O=${scratch_dir}/model_extension/stage3`,
`X=$R/code/11_Drift_vs_Selection/demographic_modeling/model_extension`, `FDA` and `SIMOUT` as above.

Stage 1, LD panel
1. `$PY $X/ld_panel/ext_mask_split.py code/11_Drift_vs_Selection/simulation $FDA/sfs/A1w.block_scan.tsv input_data/10_Introgression/peaks.bed input_data/00_Reference_Genome/${reference}.fai $W/mask`
   — the 860 callable 1-Mb blocks (`mask_blocks.tsv`) and the 13/12 calibration/held-out chromosome split
   (`chrom_split.tsv`).
2. Sample sets: `input_data/11_Drift_vs_Selection/demographic_modeling/model_extension/panel_main/samples_def.tsv`
   (name TAB comma-separated sample IDs; 14 sets: `BAM` and `TP1` first, the A1w lists in `$SIMOUT/23_inputs/`, then
   `BAM_noKin`, `TP1_noKin3` and five leave-one-out replicates of each) and `pairs_def.tsv` (BP set TAB TP set, 14
   pairs). `run_panel.sh` copies them to `$W/panel/main/`; for another run write `$W/panel/<run>/` yourself or set
   `PDEF`.
3. `bash $X/ld_panel/run_panel.sh main` (one allocation, ~120 cores) — site-table flags, per-block genotypes
   (`$W/blocks/`), the panel `panel/main/{panel.tsv,ld_curves.tsv}`, plus reproduction checks of simulation/23a, 23e,
   24a and 24b. `RUN=<name> bash $X/ld_panel/run_panel.sh variants` builds panels for other sample sets;
   `ld_panel/ext_sensitivity.py <panel.tsv ...> -o <out.tsv>` compares them with the main panel.

Stage 2b, alternative models
4. `mkdir -p $S2/fd_primary/sfs $S2/fd_primary/fits && cp $FDA/sfs/A1w.pkl $S2/fd_primary/sfs/ && cp $FDA/fits/A1w.M*.data.seed*.tsv $S2/fd_primary/fits/`
   — the drift-test spectrum and baseline fits the alternatives are fitted beside.
5. In one long job: `RUN=primary FD=$S2/fd_primary PAIRS=A1w MODELS=CB,CBT,CBTM,CBf,CBTf,CBTMf BOUNDS=primary GATE_MODELS=M0,M1,M2,M3,M0f,M1f,M2f,M3f,CB,CBT,CBTM,CBf,CBTf,CBTMf CODE=$X/stage2b IMG=$dadi_image LOGS=$S2/logs/primary bash $X/stage2b/run_stage.sh`
   — round 1, refine 2, nested pass 1, refine 3, nested pass 2, then `collect_ext.py --godambe`, `gate_ext.py`
   (`gate_ext.tsv`) and `lrt_ext.py` (`lrt_ext.tsv`).
6. `python3 $X/stage2b/write_admissible.py $S2/fd_primary/gate_ext.tsv $S2/admissible_alternatives.tsv $S2/ADMISSIBLE.DONE`
   — hand-off table for Stage 3.
7. `bash $X/stage2b/orchestrate_B.sh <wide_models|-> <var_models|->` — wide-bound refits of models with a parameter at a
   bound (`fd_wide/`) and sample-set sensitivity refits (`fd_var/`; needs the Stage 1 variant fits, see "Not included").
8. `python3 $X/stage2b/build_tables.py $S2/fd_primary --wide $S2/fd_wide --var $S2/fd_var --out <dir>` —
   `stage2b_fits.tsv`, `stage2b_lrt.tsv`, `stage2b_bound_sensitivity.tsv`, `stage2b_variant_fits.tsv`.

Stage 3, simulations (each `orch_*.sh` runs from a job work directory and submits one node job; `stage3/env.sh` sets
the paths)
9. `bash $X/stage3/orch_step1.sh $X/stage3` — `sim_index.tsv` (models, r/mu grid, seed rule), per-record site tables of
   the 860 blocks (`records/`) and their checks, the heterozygote-dropout table, the pilot reproduction and a test genome.
10. `bash $X/stage3/orch_step2.sh $X/stage3` — observed rows `primary`, `callable` and `d019` (strict-genotype mask:
    DP < 8 or GQ < 20 set missing) in `$O/obs/`.
11. `bash $X/stage3/orch_step3.sh $X/stage3 <models> <tag>` — calibration, in-sample and held-out simulation families
    over the 12-value r/mu grid, pooled into `$O/fam/fam_summary_<tag>.npz`. Baseline models (M0, M2, M3, M0f-M3f)
    first; then `bash $X/stage3/orch_alt1.sh $X/stage3` (after step 6; converts the admissible alternatives to
    `$O/alt_params.{tsv,json}` and runs the SFS unit check) and step 11 again for those alternatives.
12. `python3 $X/stage3/analysis_step3.py <fam_summary.npz> <obs_dir> $W/mask/chrom_split.tsv <models> <fits.tsv>` — r/mu
    fits per model and treatment (`obs_dir` holds the `obs_<row>_*` files that `orch_step2.sh` copies to its work dir).
13. `bash $X/stage3/orch_genomes.sh $X/stage3 <tag> 30 <model:r_true,...> 100`, then
    `python3 $X/stage3/analysis_step4.py $O/genomes/<model>/r30 <fam_summary.npz> $W/mask/chrom_split.tsv <model> <r_true> <unfolded|folded> <prefix>`
    — 100 pseudo-observed genomes per truth model and the false-rejection rates (`<prefix>_frr.tsv`).
14. `bash $X/stage3/orch_genomes.sh $X/stage3 <tag> 20 <model:r_fit,...> 100`, then
    `python3 $X/stage3/analysis_step5.py <model> $O/genomes/<model>/r20 <fam_summary.npz> <fits.tsv> <obs_dir> $W/panel/main/panel.tsv $W/mask/chrom_split.tsv <params.json> <analysis> <admissible 0|1> <p3_gating 0|1> <diag99|-> <prefix>`
    — scoring of the real data (`<prefix>_scores.tsv`; `params.json` is `simulation/params.as_run.agree.json` for
    baseline models, `$O/alt_params.json` for alternatives). `make_figures.py` holds the LD and diagnostic figure
    functions.

## Outputs
- `${scratch_dir}/model_extension/stage1/`: `mask/`, `blocks/`, `panel/<run>/{panel.tsv,ld_curves.tsv}`.
- `${scratch_dir}/model_extension/stage2b/`: `fd_primary/` (`fits/best_fits.tsv`, `gate_ext.tsv`, `lrt_ext.tsv`),
  `fd_wide/`, `fd_var/`, `admissible_alternatives.tsv`; tables from step 8.
- `${scratch_dir}/model_extension/stage3/`: `sim_index.tsv`, `records/`, `obs/`, `fam/`, `genomes/`, `alt/`; the
  r/mu fit, FRR and score tables from steps 12-14.

## Not included
The exported scripts expect these Stage 1 / pilot files, which no script in this repository writes:
- `samples_def.tsv`/`pairs_def.tsv` for panel runs other than `main`: the sets `ext_sensitivity.py` compares that
  the main run lacks (`TP1_noAdmix`, `TP1_minus10A09_r1`-`r5`, `TP1_no017`, `TP1_minus7A09_r1`-`r5`,
  `TP1_minus1BEN_r1`-`r5`; the main run's TP1 leave-one-out sets are named `TP1_minus1_r1`-`r5`). Step 2 gives the
  format.
- `$W/build2/root/output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_var2/{sfs,fits}` and
  `$W/FITS2.DONE`: spectra and baseline fits for the sample-set variants (TP1 without the admixed 2019 Apassa wobbles,
  without the 2019-017 fish, Bambomo without BAM_175), built by rerunning the drift test on those sample sets. Used only
  by the `var` track of step 7.
- `${scratch_dir}/model_extension/coord_pilot/a2/records/*.records.npz` and `a2/sim/bpp_*.npz`: coordinate-pilot
  outputs, used by `records_full.py` and `pilot_compare.py` as checks. `step1_node.sh` stops (writes
  `RECORDS_MISMATCH`) unless all 40 pilot blocks compare identical, so step 9 needs them as written.
- `$W/ld_track/observed_ld_curves.tsv`: observed LD curves read by `pilot_compare.py` (check only).
