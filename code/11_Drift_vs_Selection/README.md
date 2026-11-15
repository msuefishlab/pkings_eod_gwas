# 11_Drift_vs_Selection: drift or selection at the six BP-associated GWAS peaks

Produces: Supplementary Fig. S7 and Tables S5–S6; ESM "Tests for Selection at GWAS Peaks", "Drift Test at GWAS
Peaks" (including the exploratory regional analyses of Peak 1) and the local-ancestry follow-up in the
Supplementary Note ("Ancestry").

Each subfolder has its own README with inputs, steps and outputs. `code/`, `input_data/` and `output_data/` mirror
this layout: an analysis in `code/11_Drift_vs_Selection/<name>/` writes to `output_data/11_Drift_vs_Selection/<name>/`
and, where it has small inputs of its own, reads them from `input_data/11_Drift_vs_Selection/<name>/`.

## Subfolders

| Folder | What it does | Paper |
|---|---|---|
| `sweep_detection/` | Rsb (rehh), SweepFinder2 CLR and πBP/πTP per peak, ranked against matched non-peak regions with a two-of-three decision rule, plus the sensitivity rows | Fig. S7, Table S5 |
| `folded_clr/` | SweepFinder2 CLR on the folded spectrum (polarization check) and the six-test correction; exploratory Peak 6 test against drift-made fixed-difference clusters | Table S5 row, Fig. S7 legend |
| `demographic_modeling/` | Drift test: two-population dadi histories per BP–TP pair and E, the expected genome-wide number of SNPs as differentiated as each index SNP (`fixed_differences/`) | Table S6 |
| `demographic_modeling/model_extension/` | Bambomo histories extended with population contraction (dadi), simulated linked variation under them, and the observed LD panel they are compared with | Drift Test text (Peak 1 regional test, exploratory) |
| `simulation/` | Linked-variation tests of drift at the Bambomo peaks: simulated LD decay under the fitted histories, the empirical clustering test, the Apassa ancestry projection, and the Peaks 4/5 coalescent null and SLiM power grid they build on | Drift Test text (exploratory), Supplementary Note |
| `local_ancestry/` | Exploratory HMM local ancestry of the Bambomo-like component in the 2019 Apassa mixed fish, with simulation calibration and a peak-enrichment test | Supplementary Note, "Ancestry" |

## Order of execution
Run from the repository root after `source config.env`.

1. `sweep_detection/`: `00_fetch_inputs.sb` checks the inputs, then `ehh/00_prepare_vcfs.sb` builds the polarized
   VCFs that `folded_clr/` and `local_ancestry/` also read.
2. `demographic_modeling/fixed_differences/`: independent of step 1. Its six-species-agreement run
   (`FD_OUT=…/fixed_differences_agree`) is what `simulation/` and `folded_clr/` Part B read.
3. Then: `folded_clr/` (needs 1 and 2); `simulation/` (needs 2, and the all-sites VCF staged in 1);
   `demographic_modeling/model_extension/` (needs 2 and `simulation/`); `local_ancestry/` (needs 1 and the
   `simulation/26a`/`26c` outputs).

Large intermediates go to `${sweep_scratch}` (set in `config.env`). There is no recombination map for this
species, so any per-bp recombination rate used here is an assumption, stated in the analysis that uses it.
