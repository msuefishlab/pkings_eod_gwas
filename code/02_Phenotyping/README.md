# 02_Phenotyping: EOD waveform processing, landmark measurements and P0 calls

Produces: Fig. 1 and the per-fish EOD plots of Supplementary Data S2 (landmark analysis of averaged EODs); the
per-sweep P0 measurements and quantitative P0 calls (`p0_statistical_table.csv`) read by `14_Wobble_Analysis`.

## Inputs
- Raw EOD recordings (`input_data/02_Phenotyping/eods/`, one `.eod` or `.mat` file per fish). **Not distributed.**
  The plots made from them are Supplementary Data S2; the code is provided to document the method.
- `input_data/01_Terra/data_model/fish_data_2024a.txt` — sample sheet (columns `PN`, `EOD_FILE`, `Phenotype`);
  generate it from Supplementary Data S1 with `Rscript code/00_common/make_sample_sheets.R`.

## Steps
Run from the repository root after `source config.env`. The MATLAB steps run in one MATLAB session with
`code/02_Phenotyping` on the path; each later step uses the `normalized_eods` and `averaged_eods` structs that
`EOD_Analysis` leaves in the workspace.
1. `EOD_Analysis` (MATLAB) — reads each fish's recording (`ReadEODFileWithConversionto2019.m`), drops clipped and
   low-amplitude sweeps, baselines each sweep locally, normalizes peak-to-peak amplitude, centres on P1, and averages
   the sweeps.
2. `PlotEODs` (MATLAB; create `output_data/02_EOD_Plots/` first) — 21 landmark measurements per averaged EOD
   (`standard_eod_measurement.m`) and a four-panel PDF per fish.
3. `eod_plotter` (MATLAB) — the same four-panel plot as PNG per fish, plus `measurement_data.csv`, in
   `output_data/02_Phenotyping/`.
4. `export_per_sweep_measurements` (MATLAB) — re-measures P0 on every replicate sweep (`p0_detectors.m`) and writes
   the per-sweep P0 windows and baseline strips used by the R step.
5. `Rscript -e 'rmarkdown::render("code/02_Phenotyping/p0_statistical_detection.Rmd", output_dir = "output_data/02_Phenotyping")'`
   — label-free TP / Wobble / BP call from per-sweep P0 reproducibility, shape, dipSNR and |vP0|
   (`p0_shape_functions.R`); validated against the fish with known electrocyte anatomy. `methods.txt` describes it.
6. Optional, before step 5: `p0_phenotype_assessment.Rmd` (rendered the same way) writes `p0_confidence_table.csv`,
   which step 5 only cross-tabulates against its own call when present.

Phenotypes used in the association analyses (TP, BP, Mixed) are the manual calls in Supplementary Data S1.

## Outputs
- `output_data/02_EOD_Plots/<sample>.pdf`, `measurement_data.csv` — per-fish plots and landmarks (step 2)
- `output_data/02_Phenotyping/<sample>.png`, `measurement_data.csv` — per-fish plots and landmarks (step 3)
- `output_data/02_Phenotyping/measurement_data_per_sweep.csv`, `p0_windows*.csv`, `p0_baseline*.csv`,
  `p0_aligned_per_sweep.csv` — per-sweep measurements (step 4)
- `output_data/02_Phenotyping/p0_statistical_table.csv` and `by_stat_call/` — quantitative P0 calls (step 5)

`decodewavetex_updated.m` (B. Scott Jackson, after C. Kinnier and C. Hopkins) is a third-party helper kept with its
original credit.
