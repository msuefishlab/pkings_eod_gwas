# simulation: linked-variation tests of drift at the Bambomo peaks

Produces: the exploratory regional analyses of Peak 1 in the ESM section "Drift Test at GWAS Peaks" (simulated
LD decay under the fitted histories, 23*; the empirical clustering test, 24*), the unphased ancestry projection
of the Apassa mixed fish in the Supplementary Note (26a, 26c), and the Peaks 4/5 coalescent null and SLiM power
grid (00-07, `run_all.sh`; `METHODS.md`) whose fitted histories and workers the LD-decay scripts reuse.

## Inputs
- Drift-test outputs, `output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree/`
  (`tail/models.tsv`, `tail/index_units.tsv`, `fits/best_fits.tsv`, `checks/spectra_checks.tsv`,
  `sfs/A1w.block_scan.tsv`, `site_table/snps.*.tsv.gz`) — `demographic_modeling/fixed_differences/` run with
  `FD_OUT=…/fixed_differences_agree`. `RUN=first` reads the first-pass `fixed_differences/` run instead.
- Unphased filtered callset (`<vcf>` below),
  `output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz`
  — `03_VCF_QC_Filter`. Used by 23a, 23f, 23g, 24b, 26a and 26c.
- All-sites VCF `${sweep_scratch}/03_QC/${allsites_vcf}` — `03_VCF_QC_Filter`, placed by
  `sweep_detection/00_fetch_inputs.sb`. Stage 00 only.
- `input_data/05_PopGen/${popsfile}` (05_PopGen), `input_data/10_Introgression/peaks.bed` (10_Introgression),
  `input_data/01_Terra/data_model/fish_data_2023.txt` (`Rscript code/00_common/make_sample_sheets.R`), and the
  reference index `input_data/00_Reference_Genome/${reference}.fai` (reference genome, see the root README).
- `output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED.lmm2.assoc.txt` — 06_Association (23f).
- `params.as_run.json`, `params.as_run.agree.json` — the as-run histories, used only by
  `01_fitted_history.py --check` as a regression reference.

Software: msprime 1.4.4 / tskit 1.0.3 (pass the interpreter as `MSPRIME_PYTHON`), SLiM 5.2 (`SLIM_BIN`),
bcftools, python3 with numpy and scipy, and `images/gwas_tools.sif` for stage 00. Edit the SBATCH headers for
your cluster.

## Steps
Run from the repository root after `source config.env`. `S=code/11_Drift_vs_Selection/simulation`,
`O=output_data/11_Drift_vs_Selection/simulation`, `FD=output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree`,
and `W=${sweep_scratch}/peaks45_sim/agree` (the work directory `run_all.sh` uses).

Peaks 4/5 simulation null (`METHODS.md`):
1. `bash $S/00_extract_index_genotypes.sh` — positive control: re-extracts the Peak 1/4/5 index genotypes of
   the 114 A1w fish; every allele count must match `tail/index_units.tsv`.
2. `SLIM_BIN=<slim> MSPRIME_PYTHON=<python> bash $S/run_all.sh all` — `setup` (01: `params.json` derived
   from the drift-test fits, then `--check`), `observed` (02, 02b), `windows` (03), `coalescent` (04c; 04b for
   `RUN=first`), `slim-tasks` (05d; 05 for `RUN=first`) and `slim` (07b with `power.slim`). Stages can be run
   one at a time. 04, 05b, 05c, 06 and 07 are the first-pass variants of the workers and task generators.
   Overrides: `RUN`, `SIM_WORK`, `SIM_OUT`, `SIM_WORKERS`, `SIM_BUDGET_S`, `SIM_MSP_TASKS`, `SIM_SHARDS`,
   `SLIM_MODEL`. The summary tables `$O/01_*`–`22_*` were aggregated from these per-task outputs; no
   aggregation script is included.

LD decay under the fitted histories (needs `$W/params.json` from step 2 `setup`):
3. `mkdir -p $O/23_inputs; awk -F'\t' '$2=="BAM"{print $1}' input_data/05_PopGen/${popsfile} > $O/23_inputs/A1w_BAM.txt;
   awk -F'\t' '$2=="APA"||$2=="BEN"{print $1}' input_data/05_PopGen/${popsfile} > $O/23_inputs/A1w_TP1.txt`
   — the 67 Bambomo and 47 Apassa + Bengue fish (the same lists stage 00 writes).
4. `python3 $S/23a_observed_ld_decay.py <vcf> $O/23_inputs/A1w_<pop>.txt <pop> $FD/sfs/A1w.block_scan.tsv input_data/10_Introgression/peaks.bed $O/23a_observed_ld_<pop>.tsv`
   for `<pop>` = BAM and TP1 — observed r² decay in 200 random 1-Mb segments.
5. `mkdir -p $W/logs; sbatch -D $W --export=ALL,root=$(pwd),W=$W,MSPRIME_PYTHON=<python> -o logs/ld_decay_%A_%a.log $S/23b_sim_ld_decay.sb`
   — simulated decay, 7 histories × 12 r/μ values, into `$W/ld_sim/`.
6. `python3 $S/23e_fit_r_over_mu.py $O/23a_observed_ld_BAM.tsv $O/23a_observed_ld_TP1.tsv $W/ld_sim $O/23_r_over_mu_fit.tsv $O/23_ld_curves.tsv`
   — r/μ fit and adequacy per history.
7. Clustering controls: `python3 $S/23c_make_arm.py $W/params.json $W/cluster/ctrl_r1 1e-8` (and `ctrl_r2.5`
   with `2.5e-8`), `sbatch -D $W --export=ALL,root=$(pwd),W=$W,ARM=ctrl_r1,TASK_BASE=1,MSPRIME_PYTHON=<python> -o logs/ld_cluster_%A_%a.log $S/23c_cluster_sims.sb`
   (`TASK_BASE=101` for `ctrl_r2.5`), then `python3 $S/23d_score_clustering.py $W/cluster/ctrl_r1 $O/23d_ctrl_r1.tsv --compare $O/15_clustering_matched_window_agree.tsv 1e-08`
   (`4e-09` for `ctrl_r2.5`).
8. Diagnostics: `python3 $S/23f_dense_windows.py <vcf> $O/23f_dense_windows.tsv` (reads
   `$O/17_observed_window_counts_agree.tsv`, the stage-03 window counts); `python3 $S/23g_unlinked_floor.py <vcf> $O/23_inputs/A1w_BAM.txt $O/23a_observed_ld_BAM.tsv > $O/23g_unlinked_floor.txt`;
   `python3 $S/23h_ld_misfit_checks.py input_data/00_Reference_Genome/${reference}.fai > $O/23h_ld_misfit_checks.txt`.

Empirical clustering test (Peaks 1, 4, 5):
9. `sbatch -D $W --export=ALL,root=$(pwd),W=$W,ST=$(pwd)/$FD/site_table $S/24a_site_scan.sb` —
   `$W/24a_qualifying.tsv`, `$W/24a_density_1kb.tsv`.
10. `python3 $S/24b_empirical_clustering.py <vcf> $W/24a_qualifying.tsv $W/24a_density_1kb.tsv $O/24b_clustering`
    — block scores, matched nulls and p per peak.

Unphased ancestry projection of the Apassa fish:
11. `bash $S/26a_extract_temporal.sh <vcf> $O/26_temporal [jobs]` — DP ≥ 8, GQ ≥ 20 genotypes of Bambomo and
    Apassa fish by collection year, 2-kb regions every 20 kb outside the peaks.
12. `python3 $S/26c_ancestry_check.py <vcf> $O/26_temporal $O/26_ancestry_check.tsv` — per-fish projection
    score q onto Bambomo vs Apassa 2009 allele frequencies. It imports the loader of `26b_contemporary_ne.py`;
    26b's own contemporary-Ne estimates (`python3 $S/26b_contemporary_ne.py <vcf> $O/26_temporal $O/26_ne`)
    are not reported in the paper.

## Outputs
- `output_data/11_Drift_vs_Selection/simulation/01_*`–`22_*` — Peaks 4/5 null and power-grid tables; `${sweep_scratch}/peaks45_sim/<run>/` holds the per-task files.
- `.../23a_observed_ld_{BAM,TP1}.tsv`, `23_r_over_mu_fit.tsv`, `23_ld_curves.tsv`, `23d_ctrl_*.tsv`, `23f_dense_windows.tsv`, `23g_*.txt`, `23h_*.txt` — LD-decay test.
- `.../24b_clustering_{blocks,test}.tsv` — empirical clustering test.
- `.../26_temporal/` (genotype tables, read by `local_ancestry/` and `14_Wobble_Analysis/apassa_ancestry/`), `26_ancestry_check.tsv` — ancestry projection.
