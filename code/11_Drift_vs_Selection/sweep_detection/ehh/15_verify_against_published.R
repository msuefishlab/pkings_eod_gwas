#!/usr/bin/env Rscript
## 15_verify_against_published.R
## Check that the rehh scans in ehh/merged/ reproduce a reference scan kept in
## ehh/merged_published/ (for example the authors' scan), exactly. Optional: the pipeline skips
## this check when no reference scan is present.
##
## Why: a scan array that starts while its per-chromosome input VCFs are still being written
## reads them silently corrupted. The scan completes and the counts look right; only a
## comparison against a known-good result exposes it, so this check fails loudly.
##
## Usage (inside gwas_tools.sif, from the repo root):
##   Rscript code/11_Drift_vs_Selection/sweep_detection/ehh/15_verify_against_published.R [chr6 chr13 ...]
##
## Exit status 0 if every comparison matches exactly, 1 otherwise.

suppressPackageStartupMessages({
  library(data.table)
  library(rprojroot)
})

## .git is a directory in a clone and a file in a linked worktree, and rprojroot's
## has_file() criterion matches only the latter. Accept either, so the script runs
## wherever the repository is checked out.
root <- find_root(has_file(".git") | has_dir(".git"))
args <- commandArgs(trailingOnly = TRUE)
chroms <- if (length(args)) args else c("chr6", "chr13", "chr24")

d_new <- file.path(root, "output_data/11_Drift_vs_Selection/sweep_detection/ehh/merged")
d_pub <- file.path(root, "output_data/11_Drift_vs_Selection/sweep_detection/ehh/merged_published")
if (!dir.exists(d_pub)) stop("no published scans to compare against: ", d_pub)

read_one <- function(dir, file) {
  f <- file.path(dir, file)
  if (!file.exists(f)) return(NULL)
  x <- fread(cmd = sprintf("zcat < %s", shQuote(f)))
  setnames(x, 1:3, c("chrom", "pos", "value"))
  x[, .(chrom, pos, value)]
}

compare <- function(file, label) {
  a <- read_one(d_new, file); b <- read_one(d_pub, file)
  if (is.null(a) || is.null(b)) {
    cat(sprintf("  %-34s SKIP (missing on one side)\n", label)); return(NA)
  }
  m <- merge(a[chrom %in% chroms], b[chrom %in% chroms], by = c("chrom", "pos"),
             suffixes = c("_new", "_pub"))
  k <- m[!is.na(value_new) & !is.na(value_pub)]
  if (!nrow(k)) { cat(sprintf("  %-34s SKIP (no overlap)\n", label)); return(NA) }
  ident <- mean(k$value_new == k$value_pub)
  mx <- max(abs(k$value_new - k$value_pub))
  cat(sprintf("  %-34s n=%7d identical=%6.2f%% max|diff|=%.3g cor=%.6f\n",
              label, nrow(k), 100 * ident, mx, cor(k$value_new, k$value_pub)))
  isTRUE(all.equal(ident, 1)) && mx == 0
}

cat("Comparing the re-run scans against the published ones on: ",
    paste(chroms, collapse = ", "), "\n\n", sep = "")

results <- c(
  compare("BP1.genome_wide_ines.txt.gz", "inES BP1"),
  compare("BP2.genome_wide_ines.txt.gz", "inES BP2"),
  compare("BP3.genome_wide_ines.txt.gz", "inES BP3"),
  compare("BP1_TP1.genome_wide_rsb.txt.gz", "Rsb BP1_TP1"),
  compare("BP2_TP2.genome_wide_rsb.txt.gz", "Rsb BP2_TP2"),
  compare("BP3_TP2.genome_wide_rsb.txt.gz", "Rsb BP3_TP2")
)

tested <- results[!is.na(results)]
cat("\n")
if (!length(tested)) {
  cat("NOTHING COMPARED: no overlapping files were found.\n"); quit(status = 1)
}
if (all(tested)) {
  cat("PASS: the re-run reproduces the published scans exactly on these chromosomes.\n")
  quit(status = 0)
}
cat("FAIL: the re-run does not reproduce the published scans.\n",
    "      Do not use these scans. Check that the per-chromosome VCFs were complete and\n",
    "      unchanged before the scan array started, then re-run the scans.\n", sep = "")
quit(status = 1)
