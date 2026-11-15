# apassa_ancestry: where the Bambomo-like ancestry of the 2019 Apassa mixed fish lies (exploratory)

Produces: the genomic distribution of the Bambomo-like ancestry in the 2019 Apassa mixed ("wobble") fish
(Supplementary Note, "Ancestry"), from the phased local-ancestry fits and from phase-free projection scores.

## Inputs
- `output_data/11_Drift_vs_Selection/local_ancestry/` (`la_dir`) — local-ancestry markers, fits and simulator
  (`code/11_Drift_vs_Selection/local_ancestry`; its `la_lib.py` and `05_simulate.py` are imported)
- `output_data/11_Drift_vs_Selection/simulation/26_temporal/` (`temporal_dir`) — unphased genotypes (DP ≥ 8, GQ ≥ 20)
  from `code/11_Drift_vs_Selection/simulation/26a_extract_temporal.sh`
- `output_data/03_QC/ref_bias/ref_bias_per_sample.tsv` (`qc_tsv`) — per-sample QC from `03_VCF_QC_Filter/ref_bias_summary.R`

## Steps
Run from the repository root (python3 with numpy and scipy; about 10 minutes in total).
1. `python3 code/14_Wobble_Analysis/apassa_ancestry/01_sim_null.py <la_dir> <out>` — random-tract null: simulated fish
   at each wobble's best-match (T, m) cell plus m = 0, 50 fish per cell, seed 32032, fitted as the real fish are.
2. `python3 code/14_Wobble_Analysis/apassa_ancestry/02_analysis.py <la_dir> <temporal_dir> <qc_tsv> <out>` — Part A
   (phased): per-fish concentration and excess per chromosome and 2-Mb window over the 2009 Apassa TP baseline; Part B
   (unphased): projection scores, gates and shared-window tests.
3. `python3 code/14_Wobble_Analysis/apassa_ancestry/03_dev2.py <out> <la_dir>` — corroboration of the phased and
   unphased excess and shared-region tests against independent baselines (B9 phased, P9 unphased).

Use `out = output_data/14_Wobble_Analysis/apassa_ancestry`.

## Outputs
- `output_data/14_Wobble_Analysis/apassa_ancestry/` — `sim_null.npz` (01); `sim_null_fish.tsv`,
  `partA_per_fish.tsv`, `partA_chrom_excess.tsv`, `partB_per_fish.tsv`, `partB_chrom_excess.tsv`,
  `shared_windows.tsv`, `partB_corroboration.tsv`, `window_profiles.npz`, `results.json` (02);
  `dev2_corroboration.tsv`, `dev2_shared_windows.tsv`, `dev2_results.json` (03). The scripts write tables only; no
  figure code is included.
