# folded_clr: Peak 6 folded-spectrum CLR and the fixed-difference null

Produces: the "CLR on the folded spectrum" sensitivity row of Table S5, and the six-test Bonferroni and
Benjamini–Hochberg correction quoted in the Fig. S7 legend (`six_test_correction.tsv`, unfolded rows). Part B
(Peak 6 against drift-made fixed-difference clusters) is exploratory. `METHODS.md` describes every choice.

## Inputs
- `output_data/11_Drift_vs_Selection/sweep_detection/by_chrom/<pop>.polarized.<chr>.vcf.gz` — per-lineage polarized
  VCFs, from `sweep_detection/ehh/00_prepare_vcfs.sb` then `ehh/01b_split_by_chrom_array.sb`
- `output_data/11_Drift_vs_Selection/sweep_detection/evidence/{region_stats.tsv.gz, evidence_matrix.tsv,
  region_callable.tsv, region_differentiation.tsv, region_ines.tsv}` — from the `sweep_detection/evidence/` steps
  (02, 02b, 02c, 03)
- `output_data/11_Drift_vs_Selection/sweep_detection/sweepfinder2/sfs_input_unfolded/BP2.chr1.sf2.txt` — for the
  stage-0 check only; stage 0 regenerates it with the same converter if absent
- `output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree/{fixed_diff_cluster_counts.tsv,
  fixed_diff_clusters.*.bed}` — Part B null loci, from `demographic_modeling/fixed_differences/11_clusters.py`
- Scripts run unchanged: `sweep_detection/sweepfinder2/{vcf_to_sweepfinder2.py, 13_sf2_merge.sh}`,
  `sweep_detection/evidence/{02_region_stats.R, region_stats_lib.R}`
- Containers `images/gwas_tools.sif`, `images/sweepfinder2.sif` (see `images/README.md`)

## Steps
Run from the repository root after `source config.env`, on a compute node with singularity. `common.sh` sets
every path; `bash code/11_Drift_vs_Selection/folded_clr/run_all.sh` runs stages 0–8 in order.

| Stage | Script | Produces |
|---|---|---|
| 0 | `00_folded_precondition.sh` (`check_folded_conversion.py`) | `audit/folded_conversion_check.tsv`; blocking (METHODS A3) |
| 1 | `01_convert.sh` → `01_convert_worker.sh` | `sweepfinder2/sfs_input_folded/<pop>.<chr>.sf2.txt` (125 files) |
| 2 | `02_background.sh` → `02_background_worker.sh` | `sweepfinder2/background_folded/<pop>.spectrum` |
| 3 | `03_scan.sh` → `03_scan_worker.sh` | `sweepfinder2/clr_folded/<pop>.<chr>.clr.txt`, `folded_scan_gridpoints.tsv` |
| 4 | `04_merge_folded.sh` (calls `sweep_detection/sweepfinder2/13_sf2_merge.sh folded`) | `sweepfinder2/merged/<pop>.genome_wide_clr_folded.txt` |
| 5 | `05_region_stats.sh` | `evidence/region_stats_folded.tsv.gz`, `evidence/region_stats_operators_folded.tsv` |
| 6 | `06_evidence_matrix.sh` → `06_evidence_matrix_clrstat.R` | `evidence_matrix_{rederived_unfolded,folded}.tsv`, `calibration_*.tsv` |
| 7 | `07_phaseB.sh` → `07_phaseB.R` | `phaseB_precondition.tsv`, `phaseB_percentiles.tsv`, `phaseB_cluster_values.tsv.gz` |
| 8 | `08_correction.sh` → `08_correction.R` | `reproduction_check.tsv` (blocking), `folded_vs_unfolded_percentiles.tsv`, `six_test_correction.tsv`, `peak6_folded_verdict.tsv` |

`sweep_detection/sweepfinder2/10_sf2_convert.sb` with `MODE=folded` writes the same folded inputs as stage 1
(same converter, same flags).

Preconditions: stage 0 (folded conversion is polarization-free) and stage 8 (the re-derived unfolded run
reproduces the headline evidence matrix to 1e-6) stop the analysis if they fail. In Part B, Peak 1 / BP1 (pair
A1w) has fewer than 200 drift-made clusters and is recorded as a precondition failure in
`phaseB_precondition.tsv`; Part B is reported for Peak 6 only.

## Outputs
- `output_data/11_Drift_vs_Selection/folded_clr/` — the tables listed above, `audit/`, `logs/`
- `output_data/11_Drift_vs_Selection/sweep_detection/sweepfinder2/{sfs_input_folded,background_folded,clr_folded,merged}/`
  — SweepFinder2 intermediates, where `region_stats_lib.R` reads them
