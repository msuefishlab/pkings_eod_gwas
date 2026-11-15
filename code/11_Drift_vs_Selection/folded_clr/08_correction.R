#!/usr/bin/env Rscript
## 08_correction.R
## Folded vs unfolded side by side for all six peaks, and the six-test correction.
##
## The self-check below is blocking: the re-derived unfolded percentiles must reproduce the
## headline evidence matrix (sweep_detection/evidence/evidence_matrix.tsv) to 1e-6, or nothing
## downstream is reported.
suppressPackageStartupMessages({ library(data.table); library(rprojroot) })
root <- find_root(has_file(".git") | has_dir(".git"))
args <- commandArgs(trailingOnly = TRUE)
arg_of <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
OUTD <- arg_of("--out_dir", Sys.getenv("FCLR_OUT",
                file.path(root, "output_data/11_Drift_vs_Selection/folded_clr")))
if (!startsWith(OUTD, "/")) OUTD <- file.path(root, OUTD)
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
dir_ev <- file.path(root, "output_data/11_Drift_vs_Selection/sweep_detection/evidence")

## Re-derived unfolded and folded runs are both written by 06_evidence_matrix_clrstat.R into
## OUTD; the reference ("published") matrix is the headline one in the sweep-detection evidence directory.
em_u <- fread(arg_of("--rederived", file.path(OUTD, "evidence_matrix_rederived_unfolded.tsv")))
em_f <- fread(arg_of("--folded",    file.path(OUTD, "evidence_matrix_folded.tsv")))
pub  <- fread(arg_of("--published", file.path(dir_ev, "evidence_matrix.tsv")))

## -------- self-check: the re-derived unfolded run must reproduce the headline numbers
chk <- merge(em_u[statistic %in% c("rsb","clr","pi_ratio"), .(peak, statistic, peak_value, pct_matched, empirical_p)],
             pub[statistic %in% c("rsb","clr","pi_ratio"), .(peak, statistic, peak_value, pct_matched, empirical_p)],
             by = c("peak","statistic"), suffixes = c("_rederived","_published"))
chk[, `:=`(d_value = abs(peak_value_rederived - peak_value_published),
           d_pct   = abs(pct_matched_rederived - pct_matched_published),
           d_p     = abs(empirical_p_rederived - empirical_p_published))]
chk[, reproduced := d_value < 1e-6 & d_pct < 1e-6 & d_p < 1e-6]
fwrite(chk, file.path(OUTD, "reproduction_check.tsv"), sep = "\t")
cat("--- reproduction check: re-derived unfolded vs committed evidence matrix ---\n")
cat(sprintf("rows %d, reproduced %d, max |d pct_matched| %.3e, max |d p| %.3e\n",
            nrow(chk), sum(chk$reproduced), max(chk$d_pct), max(chk$d_p)))
if (!all(chk$reproduced)) { print(chk[reproduced == FALSE]); stop("the re-derived unfolded run does not reproduce the committed evidence matrix") }

## -------- per-peak, per-statistic percentile table
long <- rbind(
  em_u[, .(peak, lineages, statistic, direction, peak_value, pct_matched, pct_unmatched,
           tp_value, tp_pct, n_matched_null, run = "unfolded_primary")],
  em_f[, .(peak, lineages, statistic, direction, peak_value, pct_matched, pct_unmatched,
           tp_value, tp_pct, n_matched_null, run = "folded_primary")])
stat_tbl <- dcast(unique(long[, .(peak, lineages, statistic, peak_value, pct_matched)]),
                  peak + lineages ~ statistic, value.var = c("peak_value", "pct_matched"))
fwrite(stat_tbl, file.path(OUTD, "folded_vs_unfolded_percentiles.tsv"), sep = "\t")
cat("\n--- folded vs unfolded percentiles ---\n"); print(stat_tbl)

## -------- six-test correction, both runs
calls <- function(em, lab, clr_stat) {
  d <- unique(em[, .(peak, lineages, composite, supported, strong, empirical_p)])
  setorder(d, peak)
  d[, `:=`(run = lab, clr_statistic = clr_stat,
           bonferroni = pmin(1, 6 * empirical_p),
           bh_q = p.adjust(empirical_p, method = "BH"))]
  d[, `:=`(survives_bonferroni_05 = bonferroni <= 0.05, survives_bh_05 = bh_q <= 0.05)]
  d[]
}
corr <- rbind(calls(em_u, "unfolded_primary", "clr"),
              calls(em_f, "folded_primary", "clr_folded"))
setorder(corr, peak, run)
fwrite(corr, file.path(OUTD, "six_test_correction.tsv"), sep = "\t")
cat("\n--- six-test correction ---\n"); print(corr)

p6 <- corr[peak == 6]
f <- p6[run == "folded_primary"]
verdict <- if (!f$supported) "NOT ROBUST" else if (f$bh_q <= 0.05) "ROBUST" else "WEAKENED"
cat(sprintf("\nPeak 6 folded verdict (rule fixed in METHODS.md A6): %s\n", verdict))
fwrite(data.table(peak = 6, verdict = verdict, supported_folded = f$supported,
                  strong_folded = f$strong, composite_folded = f$composite,
                  empirical_p_folded = f$empirical_p, bonferroni_folded = f$bonferroni,
                  bh_q_folded = f$bh_q,
                  composite_unfolded = p6[run == "unfolded_primary", composite],
                  empirical_p_unfolded = p6[run == "unfolded_primary", empirical_p],
                  bh_q_unfolded = p6[run == "unfolded_primary", bh_q]),
       file.path(OUTD, "peak6_folded_verdict.tsv"), sep = "\t")
