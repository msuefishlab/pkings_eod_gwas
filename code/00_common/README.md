# 00_common: shared helpers

Produces: the sample sheets every module reads, and the R Markdown rendering wrappers other modules call.

## Inputs
- `input_data/01_Terra/data_model/SupplementalDataS1_Fish_Specimens.xlsx` — Supplementary Data S1, downloaded from the
  paper

## Steps
Run from the repository root after `source config.env`.
1. `Rscript code/00_common/make_sample_sheets.R` — builds the sample sheets from Supplementary Data S1. Run this first.

The other scripts are called by module scripts, not run on their own:
- `render_rmd_report.R DATA_PATH OUTPUT_PREFIX REPORT.Rmd [MAF_THRESHOLD]` — renders an R Markdown report that
  declares the params `data_path` and `MAF_thresh` (used by `03_VCF_QC_Filter`, `04_Phylogenomics`, `08_Peak_Analysis`).
- `render_gwas_report.R DATA_PATH OUTPUT_PREFIX REPORT.Rmd MAF_THRESHOLD [TEST]` — the same, plus the param `test`
  (`lrt`, the default and the published run, or `wald`) (used by `06_Association` and `08_Peak_Analysis`).

## Outputs
- `input_data/01_Terra/data_model/fish_data_2023.txt`, `fish_data_2024a.txt` — sample sheets (git ignores them)
- `OUTPUT_PREFIX.report.html` — from the render wrappers
