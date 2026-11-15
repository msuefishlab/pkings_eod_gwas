#!/usr/bin/env Rscript
# ============================================================================
# 07_make_p0_phenos.R
#
# Phenotype files for the Bambomo (BAM) wobble scans with the label-free P0
# classifier as well as the manual Wobble/BP labels.
#
# The classifier (code/02_Phenotyping/p0_statistical_detection.Rmd) re-calls the
# 25 manual BAM wobbles as 15 Wobble / 4 BP / 5 Ambiguous / 1 TP, and the
# anatomy-confirmed wobble BAM_6494 carries the manual label BP. This script writes one
# multi-column GEMMA phenotype file (select a column with `gemma -n`) aligned to
# the existing BAM genotype file, so the geno, map and kinship of
# output_data/06_Association/BAM_ONLY_BP0_TPWOB1/ (phased VCF, MAF 0.05 on the
# 67 fish) are reused unchanged.
#
# Columns (NA = excluded from that scan):
#   1 MANUAL_BIN        manual labels: BP=0, Wobble=1 (TP NA). Baseline.
#   2 CLS_BIN           classifier call: BP=0, Wobble=1; Ambiguous/TP NA.
#   3 CLS_STRICT        manual and classifier agree: BP/BP=0, Wobble/Wobble=1.
#   4 REPRO_INT         rank-INT of `repro` (P0 reproducibility; NA if <4 sweeps).
#   5 VP0_INT           rank-INT of |vP0| (P0 depth relative to the whole EOD).
#   6 REPRO_INT_FULLBP  as 4, only fish carrying a BP haplotype at peaks 1, 4 and 5.
#   7 VP0_INT_FULLBP    as 5, same restriction (the modifier scan).
# BAM_6494 is anatomy-confirmed Wobble: coded 1 in columns 2-3 (classifier says
# Ambiguous; its manual label BP is kept in column 1).
# The one BAM TP fish (BAM_6597) is NA everywhere.
#
# Inputs (relative to repo root):
#   output_data/06_Association/BAM_ONLY_BP0_TPWOB1/BAM_ONLY_BP0_TPWOB1.order
#     (sample order of the GEMMA files; 06_Association/02_run_rel_mat.sh)
#   output_data/02_Phenotyping/p0_statistical_table.csv
#   output_data/14_Wobble_Analysis/PKINGS_WOB_MULTILOCUS/multilocus_peak_haplotype_wide.csv
#   input_data/01_Terra/data_model/fish_data_2023.txt
#
# Outputs (input_data/14_Wobble_Analysis/BAM_P0/; shipped with the repository, because
# p0_statistical_table.csv needs the raw EOD recordings):
#   BAM_P0.pheno           7 columns, no header, rows in .order order
#   BAM_P0.pheno_columns.tsv  column index, name, n, cases/controls
#   BAM_P0.samples.tsv     per-fish labels and traits (shared table for later steps)
#
# Run:  Rscript code/14_Wobble_Analysis/07_make_p0_phenos.R
# ============================================================================

suppressPackageStartupMessages(library(tidyverse))

root <- rprojroot::find_root(rprojroot::is_git_root)
setwd(root)

outdir <- file.path("input_data", "14_Wobble_Analysis", "BAM_P0")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ANATOMY_WOBBLE <- "BAM_6494"

order_ids <- readLines("output_data/06_Association/BAM_ONLY_BP0_TPWOB1/BAM_ONLY_BP0_TPWOB1.order")

meta <- read_tsv("input_data/01_Terra/data_model/fish_data_2023.txt", show_col_types = FALSE) %>%
  select(PN = `entity:participant_id`, manual = Phenotype, SL, Sex)
p0 <- read_csv("output_data/02_Phenotyping/p0_statistical_table.csv", show_col_types = FALSE) %>%
  select(PN = sample_name, n_sweeps, repro, vP0, p_present, presence_tier, stat_call, confidence)
hap <- read_csv("output_data/14_Wobble_Analysis/PKINGS_WOB_MULTILOCUS/multilocus_peak_haplotype_wide.csv",
                show_col_types = FALSE) %>%
  select(PN, bpN_1, bpN_4, bpN_5)

d <- tibble(PN = order_ids) %>%
  left_join(meta, by = "PN") %>%
  left_join(p0, by = "PN") %>%
  left_join(hap, by = "PN")

stopifnot(nrow(d) == length(order_ids), !anyNA(d$manual), !anyNA(d$stat_call), !anyNA(d$bpN_1))

# Rank-based inverse normal transform over the non-NA values.
rint <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  out[ok] <- qnorm((rank(x[ok], ties.method = "average") - 0.5) / sum(ok))
  out
}

d <- d %>%
  mutate(
    abs_vP0  = abs(vP0),
    is_TP    = manual == "TP",
    full_bp  = bpN_1 >= 1 & bpN_4 >= 1 & bpN_5 >= 1,
    MANUAL_BIN = case_when(manual == "BP" ~ 0, manual == "Wobble" ~ 1),
    CLS_BIN    = case_when(PN == ANATOMY_WOBBLE ~ 1,
                           stat_call == "BP" ~ 0, stat_call == "Wobble" ~ 1),
    CLS_STRICT = case_when(PN == ANATOMY_WOBBLE ~ 1,
                           manual == "BP" & stat_call == "BP" ~ 0,
                           manual == "Wobble" & stat_call == "Wobble" ~ 1),
    REPRO_INT  = rint(if_else(is_TP, NA_real_, repro)),
    VP0_INT    = rint(if_else(is_TP, NA_real_, abs_vP0)),
    REPRO_INT_FULLBP = rint(if_else(is_TP | !full_bp, NA_real_, repro)),
    VP0_INT_FULLBP   = rint(if_else(is_TP | !full_bp, NA_real_, abs_vP0))
  )

pheno_cols <- c("MANUAL_BIN", "CLS_BIN", "CLS_STRICT", "REPRO_INT", "VP0_INT",
                "REPRO_INT_FULLBP", "VP0_INT_FULLBP")

d %>%
  select(all_of(pheno_cols)) %>%
  mutate(across(everything(), ~ if_else(is.na(.x), "NA", format(.x, digits = 6, trim = TRUE)))) %>%
  write_delim(file.path(outdir, "BAM_P0.pheno"), delim = "\t", col_names = FALSE)

col_key <- map_dfr(seq_along(pheno_cols), function(i) {
  x <- d[[pheno_cols[i]]]
  binary <- all(x[!is.na(x)] %in% c(0, 1))
  tibble(col = i, name = pheno_cols[i], n = sum(!is.na(x)),
         n_case = if (binary) sum(x == 1, na.rm = TRUE) else NA_integer_,
         n_ctrl = if (binary) sum(x == 0, na.rm = TRUE) else NA_integer_)
})
write_tsv(col_key, file.path(outdir, "BAM_P0.pheno_columns.tsv"))

d %>%
  select(PN, manual, stat_call, presence_tier, confidence, n_sweeps, repro, vP0, abs_vP0,
         p_present, bpN_1, bpN_4, bpN_5, full_bp, SL, Sex, all_of(pheno_cols)) %>%
  write_tsv(file.path(outdir, "BAM_P0.samples.tsv"))

print(col_key, n = Inf)
