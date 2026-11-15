# demographic_modeling: dadi histories and the drift test at the GWAS index SNPs

Produces: Table S6 and the ESM Methods section "Drift Test at GWAS Peaks"; the exploratory regional test of Peak 1
under extended histories in the same section.

- `fixed_differences/` — the drift test: two-population dadi fits per BP–TP pair, the expected and observed number
  of SNPs as differentiated as each index SNP, and Table S6 (`13_report_tables.py`). Run this first; its README
  gives the steps and `METHODS.md` the definitions.
- `model_extension/` — exploratory: alternative histories with population contraction fitted to the Bambomo–TP1
  pair, and simulations of linked variation under the fitted histories to compare linkage-disequilibrium decay
  with the observed panel. Reads the outputs of `fixed_differences/`; see its README.

Inputs and outputs are under `input_data/11_Drift_vs_Selection/demographic_modeling/` and
`output_data/11_Drift_vs_Selection/demographic_modeling/`.
