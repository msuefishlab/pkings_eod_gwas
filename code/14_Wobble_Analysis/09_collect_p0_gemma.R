#!/usr/bin/env Rscript
# ============================================================================
# 09_collect_p0_gemma.R
#
# Summarizes the BAM_P0 scans (08_run_p0_gemma.sh) from the reduced files made
# on the cluster by 09_reduce_p0_gemma.sh (<NAME>.top.tsv = SNPs with p < 1e-3,
# <NAME>.perm_min.tsv = min-p line per permutation), copied to
# output_data/14_Wobble_Analysis/BAM_P0/<NAME>/.
#
# Per scan: permutation thresholds (95th = significant, 67th = suggestive, the
# pipeline convention), peaks (suggestive SNPs merged within 75 kb), the best
# p inside main peaks 1/4/5, and p at the recurrent wobble leads.
#
# Outputs (output_data/14_Wobble_Analysis/BAM_P0/SUMMARY/):
#   scan_summary.tsv   one row per scan: n, thresholds, top SNP, n peaks
#   scan_peaks.tsv     suggestive/significant peaks per scan
#   scan_landmarks.tsv best p within main peaks 1/4/5 and at recurrent leads
#   scan_manhattan.png p < 1e-3 SNPs per scan with thresholds
#
# Run:  Rscript code/14_Wobble_Analysis/09_collect_p0_gemma.R
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

root <- rprojroot::find_root(rprojroot::is_git_root)
setwd(root)

base   <- file.path("output_data", "14_Wobble_Analysis", "BAM_P0")
outdir <- file.path(base, "SUMMARY")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cols <- read_tsv("input_data/14_Wobble_Analysis/BAM_P0/BAM_P0.pheno_columns.tsv", show_col_types = FALSE)

# Main peaks of PKINGS_ALL_WOB_EXCLUDED with Bambomo origin (Peaks 1, 4, 5).
main_peaks <- read_tsv("output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed",
                       col_names = c("chrnum", "start", "end", "peak", "logp", "index_pos"),
                       show_col_types = FALSE) %>%
  filter(peak %in% c(1, 4, 5)) %>%
  transmute(landmark = paste0("main_peak_", peak), chr = paste0("chr", chrnum), start, end)

# Recurrent leads from earlier wobble scans (10_recurrent_leads.Rmd), +-25 kb.
leads <- tribble(
  ~landmark,          ~chr,    ~pos,
  "lead_chr9_31.19",  "chr9",  31189837,
  "lead_chr4_19.24",  "chr4",  19235836,
  "lead_chr13_31.60", "chr13", 31603593,
  "lead_chr3_17.49",  "chr3",  17486697,
  "lead_chr7_13.39",  "chr7",  13388648
) %>% transmute(landmark, chr, start = pos - 25000, end = pos + 25000)
landmarks <- bind_rows(main_peaks, leads)

merge_peaks <- function(hits, gap = 75000) {
  if (nrow(hits) == 0) {
    return(tibble(chr = character(), start = numeric(), end = numeric(), n_snps = integer(),
                  index_snp = character(), index_ps = numeric(), min_p = numeric(), index_af = numeric()))
  }
  hits %>%
    arrange(chr, ps) %>%
    group_by(chr) %>%
    mutate(block = cumsum(c(1, diff(ps) > gap))) %>%
    group_by(chr, block) %>%
    summarise(start = min(ps), end = max(ps), n_snps = n(),
              index_snp = rs[which.min(p_wald)], index_ps = ps[which.min(p_wald)],
              min_p = min(p_wald), index_af = af[which.min(p_wald)], .groups = "drop") %>%
    select(-block)
}

summ <- list(); peaks <- list(); marks <- list(); man <- list()
for (nm in cols$name) {
  top_f  <- file.path(base, nm, paste0(nm, ".top.tsv"))
  perm_f <- file.path(base, nm, paste0(nm, ".perm_min.tsv"))
  if (!file.exists(top_f)) { message("missing ", top_f); next }
  top <- fread(top_f) %>% as_tibble()
  perm <- if (file.exists(perm_f)) fread(perm_f) %>% as_tibble() else tibble(p_wald = numeric())
  perm_lp <- -log10(perm$p_wald)
  sig  <- if (length(perm_lp) >= 20) unname(quantile(perm_lp, 0.95)) else NA_real_
  sugg <- if (length(perm_lp) >= 20) unname(quantile(perm_lp, 0.67)) else NA_real_

  pk <- merge_peaks(filter(top, -log10(p_wald) >= sugg)) %>%
    mutate(scan = nm, level = if_else(-log10(min_p) >= sig, "significant", "suggestive"))
  peaks[[nm]] <- pk

  marks[[nm]] <- landmarks %>%
    rowwise() %>%
    mutate(scan = nm,
           best_p = { x <- top$p_wald[top$chr == chr & top$ps >= start & top$ps <= end]
                      if (length(x)) min(x) else NA_real_ }) %>%
    ungroup() %>%
    mutate(best_p_label = if_else(is.na(best_p), ">1e-3", formatC(best_p, format = "e", digits = 2)))

  best <- top %>% slice_min(p_wald, n = 1, with_ties = FALSE)
  summ[[nm]] <- tibble(scan = nm, n = cols$n[cols$name == nm],
                       n_perm = length(perm_lp), thr_sig = sig, thr_sugg = sugg,
                       top_snp = best$rs, top_logp = -log10(best$p_wald), top_af = best$af,
                       n_snps_p1e3 = nrow(top),
                       n_sig_peaks = sum(pk$level == "significant"),
                       n_sugg_peaks = nrow(pk))
  man[[nm]] <- top %>% transmute(scan = nm, chr, ps, logp = -log10(p_wald), thr_sig = sig, thr_sugg = sugg)
}

summ_tbl <- bind_rows(summ)
write_tsv(summ_tbl, file.path(outdir, "scan_summary.tsv"))
write_tsv(bind_rows(peaks), file.path(outdir, "scan_peaks.tsv"))
marks_tbl <- bind_rows(marks)
write_tsv(marks_tbl, file.path(outdir, "scan_landmarks.tsv"))

chr_levels <- paste0("chr", 1:25)
scan_lab <- summ_tbl %>%
  mutate(lab = sprintf("%s\nn = %d, sig %.1f", scan, as.integer(n), thr_sig)) %>%
  select(scan, lab)
man_tbl <- bind_rows(man) %>%
  left_join(scan_lab, by = "scan") %>%
  mutate(chr = factor(chr, levels = chr_levels),
         lab = factor(lab, levels = scan_lab$lab[match(cols$name, scan_lab$scan)]),
         shade = factor(as.integer(chr) %% 2))
thr_long <- man_tbl %>%
  distinct(lab, thr_sig, thr_sugg) %>%
  pivot_longer(c(thr_sig, thr_sugg), names_to = "threshold", values_to = "y") %>%
  mutate(threshold = recode(threshold, thr_sig = "significant (95th pct)", thr_sugg = "suggestive (67th pct)"))
p <- ggplot(man_tbl, aes(ps / 1e6, logp)) +
  geom_point(aes(colour = shade), size = 0.6, alpha = 0.5, shape = 16, stroke = 0, show.legend = FALSE) +
  geom_hline(data = thr_long, aes(yintercept = y, linetype = threshold), colour = "firebrick", linewidth = 0.5) +
  scale_colour_manual(values = c("0" = "grey20", "1" = "steelblue4")) +
  scale_linetype_manual(values = c("significant (95th pct)" = "solid", "suggestive (67th pct)" = "dashed"),
                        name = "Permutation threshold") +
  facet_grid(lab ~ chr, scales = "free_x", space = "free_x",
             labeller = labeller(chr = function(x) sub("chr", "", x))) +
  labs(x = "Chromosome", y = expression(-log[10](p[Wald])),
       title = "BAM wobble scans: manual vs classifier labels and continuous P0 traits (SNPs with p < 1e-3)") +
  theme_bw(base_size = 12) +
  theme(panel.spacing.x = unit(1, "pt"), panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 9), strip.text.y = element_text(angle = 0, size = 10),
        legend.position = "bottom")
ggsave(file.path(outdir, "scan_manhattan.png"), p, width = 18, height = 12, dpi = 150)

print(summ_tbl, width = Inf)
marks_tbl %>% select(scan, landmark, best_p_label) %>%
  pivot_wider(names_from = landmark, values_from = best_p_label) %>% print(width = Inf)
