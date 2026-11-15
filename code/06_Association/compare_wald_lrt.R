#!/usr/bin/env Rscript
# Wald vs LRT comparison for the main GWAS (PKINGS_ALL_WOB_EXCLUDED; ESM Association Analysis).
# Compares three runs: Wald (-lmm 1, phased VCF), a matched LRT (-lmm 2) on the same phased inputs and the
# same 100 permuted phenotypes, and the main LRT run (unphased VCF; PEAKS.bed + min.txt only).
# Wald vs matched LRT isolates the test; matched LRT vs published LRT isolates the VCF change.
#
# Rscript compare_wald_lrt.R OUTDIR
# OUTDIR must hold published_PEAKS.bed, published_min.txt, wald_perm_min.tsv and lrt_perm_min.tsv
# (perm_id <tab> min p, written by compare_wald_lrt.sb).

suppressMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
root <- rprojroot::find_root(".git/index")
assoc_dir <- file.path(root, "output_data", "06_Association")
MAF <- 0.15

runs <- list(
  wald = list(assoc = file.path(assoc_dir, "PKINGS_ALL_WOB_EXCLUDED", "PKINGS_ALL_WOB_EXCLUDED.lmm1.assoc.txt"),
              min = file.path(assoc_dir, "PKINGS_ALL_WOB_EXCLUDED", "PKINGS_ALL_WOB_EXCLUDED_permution", "min.txt"),
              peaks = file.path(assoc_dir, "PKINGS_ALL_WOB_EXCLUDED", "PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed"),
              pcol = "p_wald", mincol = "V12"),
  lrt_phased = list(assoc = file.path(assoc_dir, "PKINGS_ALL_WOB_EXCLUDED_LRT", "PKINGS_ALL_WOB_EXCLUDED_LRT.lmm2.assoc.txt"),
                    min = file.path(assoc_dir, "PKINGS_ALL_WOB_EXCLUDED_LRT", "PKINGS_ALL_WOB_EXCLUDED_LRT_permution", "min.txt"),
                    peaks = file.path(assoc_dir, "PKINGS_ALL_WOB_EXCLUDED_LRT", "PKINGS_ALL_WOB_EXCLUDED_LRT_PEAKS.bed"),
                    pcol = "p_lrt", mincol = "V10"),
  lrt_published = list(min = file.path(outdir, "published_min.txt"),
                       peaks = file.path(outdir, "published_PEAKS.bed"),
                       mincol = "V10")
)

# same rule as gwas_report.Rmd: 95th / 67th percentile of the permutation min p
thresholds <- function(min_file, col) {
  p <- sort(fread(min_file, header = FALSE)[[col]])
  n <- length(p)
  c(n_perm = n, sig = -log10(p[n - floor(n * 0.95)]), sugg = -log10(p[n - floor(n * 0.67)]))
}

read_peaks <- function(f) {
  b <- fread(f, header = FALSE, col.names = c("chr", "start", "end", "peak", "max_log_p", "top_bp"))
  GRanges(paste0("chr", b$chr), IRanges(b$start, b$end), peak = b$peak, max_log_p = b$max_log_p, top_bp = b$top_bp)
}

read_assoc <- function(f, pcol) {
  d <- fread(f, select = c("chr", "ps", "af", pcol))
  setnames(d, pcol, "p")
  d[, log_p := -log10(p)]
  d
}

thr <- lapply(runs, function(r) thresholds(r$min, r$mincol))
peaks <- lapply(runs, function(r) read_peaks(r$peaks))
pub <- peaks$lrt_published

assoc <- lapply(runs[c("wald", "lrt_phased")], function(r) read_assoc(r$assoc, r$pcol))

lambda_gc <- function(p) qchisq(median(p, na.rm = TRUE), 1, lower.tail = FALSE) / qchisq(0.5, 1)

summary_tab <- rbindlist(lapply(names(runs), function(k) {
  a <- assoc[[k]]
  data.table(run = k,
             n_snps_tested = if (is.null(a)) NA_integer_ else nrow(a),
             n_snps_af_filter = if (is.null(a)) NA_integer_ else a[af > MAF, .N],
             lambda_gc = if (is.null(a)) NA_real_ else round(lambda_gc(a$p), 3),
             n_perm = thr[[k]][["n_perm"]],
             sig = round(thr[[k]][["sig"]], 2), sugg = round(thr[[k]][["sugg"]], 2),
             n_peaks = length(peaks[[k]]))
}))
fwrite(summary_tab, file.path(outdir, "summary.tsv"), sep = "\t")

# each published peak, scored in each new run: max -log10 p inside the published interval, SNPs above that
# run's suggestive threshold, and whether that run's own peak calling overlaps it. "af" uses the pipeline's
# af > 0.15 filter; "maf" uses the symmetric min(af, 1 - af) > 0.15 that the Methods describe.
score_region <- function(a, gr, sugg) {
  rbindlist(lapply(seq_along(gr), function(i) {
    w <- a[chr == as.character(seqnames(gr))[i] & ps >= start(gr)[i] & ps <= end(gr)[i]]
    wa <- w[af > MAF]
    wm <- w[pmin(af, 1 - af) > MAF]
    data.table(max_log_p_af = round(max(c(wa$log_p, 0)), 2), n_above_sugg_af = sum(wa$log_p > sugg),
               max_log_p_maf = round(max(c(wm$log_p, 0)), 2), n_above_sugg_maf = sum(wm$log_p > sugg))
  }))
}

pub_check <- data.table(peak = paste0("P", pub$peak), chr = as.character(seqnames(pub)),
                        start = start(pub), end = end(pub), published_max_log_p = round(pub$max_log_p, 2))
for (k in c("wald", "lrt_phased")) {
  s <- score_region(assoc[[k]], pub, thr[[k]][["sugg"]])
  setnames(s, paste0(k, "_", names(s)))
  hit <- findOverlaps(pub, peaks[[k]])
  called <- rep(NA_character_, length(pub))
  called[queryHits(hit)] <- paste0("peak ", peaks[[k]]$peak[subjectHits(hit)])
  pub_check <- cbind(pub_check, s, setNames(data.table(called), paste0(k, "_called")))
}
fwrite(pub_check, file.path(outdir, "published_peaks_check.tsv"), sep = "\t")

new_peaks <- rbindlist(lapply(c("wald", "lrt_phased"), function(k) {
  g <- peaks[[k]][!overlapsAny(peaks[[k]], pub)]
  if (length(g) == 0) return(NULL)
  data.table(run = k, chr = as.character(seqnames(g)), start = start(g), end = end(g),
             max_log_p = round(g$max_log_p, 2), top_bp = g$top_bp)
}))
fwrite(new_peaks, file.path(outdir, "new_peaks_not_in_published.tsv"), sep = "\t")

# paired permutation minima: the same permuted phenotype under each test
perm <- merge(fread(file.path(outdir, "wald_perm_min.tsv"), col.names = c("perm", "wald_min_p")),
              fread(file.path(outdir, "lrt_perm_min.tsv"), col.names = c("perm", "lrt_min_p")), by = "perm")
perm[, `:=`(wald = -log10(wald_min_p), lrt = -log10(lrt_min_p))]
fwrite(perm, file.path(outdir, "perm_min_paired.tsv"), sep = "\t")

# per-SNP agreement between the two tests on identical SNPs
m <- merge(assoc$wald[, .(chr, ps, af, wald = log_p)], assoc$lrt_phased[, .(chr, ps, lrt = log_p)], by = c("chr", "ps"))
snp_cor <- m[af > MAF, .(n = .N, spearman = cor(wald, lrt, method = "spearman"),
                         median_diff_wald_minus_lrt = median(wald - lrt))]
top <- m[af > MAF & (wald > 4 | lrt > 4)]
snp_cor[, `:=`(n_top = nrow(top), median_diff_top = median(top$wald - top$lrt))]
fwrite(snp_cor, file.path(outdir, "snp_agreement.tsv"), sep = "\t")

p1 <- ggplot(top, aes(lrt, wald)) + geom_point(size = 0.4, alpha = 0.3) +
  geom_abline(linetype = 2, colour = "grey50") +
  geom_hline(yintercept = thr$wald[["sugg"]], colour = "steelblue") +
  geom_vline(xintercept = thr$lrt_phased[["sugg"]], colour = "steelblue") +
  labs(x = "LRT -log10 p", y = "Wald -log10 p", title = "SNPs with -log10 p > 4 in either test",
       subtitle = "blue: each test's suggestive threshold") + theme_bw(base_size = 12)
p2 <- ggplot(perm, aes(lrt, wald)) + geom_point() + geom_abline(linetype = 2, colour = "grey50") +
  labs(x = "LRT min -log10 p", y = "Wald min -log10 p", title = "Genome-wide minimum per permutation",
       subtitle = paste0("n = ", nrow(perm), " paired permuted phenotypes")) + theme_bw(base_size = 12)
ggsave(file.path(outdir, "wald_vs_lrt_snps.png"), p1, width = 6, height = 5, dpi = 200)
ggsave(file.path(outdir, "wald_vs_lrt_perm_min.png"), p2, width = 6, height = 5, dpi = 200)

print(summary_tab)
print(pub_check)
print(new_peaks)
print(snp_cor)
