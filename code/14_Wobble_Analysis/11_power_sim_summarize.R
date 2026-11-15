#!/usr/bin/env Rscript
# ============================================================================
# 11_power_sim_summarize.R
#
# Scores the power simulation (11_power_sim_make.R -> power_sim.sb) against the
# genome-wide permutation thresholds of the real BAM_P0 scans
# (95th percentile of permutation minimum p, 09_collect_p0_gemma.R):
#   CLS -> CLS_BIN; MANUAL, MANUAL_CLEAN -> MANUAL_BIN; QUANT -> VP0_INT.
#
# Inputs (output_data/14_Wobble_Analysis/BAM_P0/):
#   POWER_SIM/sim_design.tsv, POWER_SIM/results/chunk_*.tsv (from power_sim.sb)
#   SUMMARY/scan_summary.tsv
# Outputs (POWER_SIM/):
#   power_table.tsv   power per design x MAF bin x effect (+ nominal p<0.05 rate)
#   power_80.tsv      smallest simulated effect with power >= 0.8 per design x MAF bin
#   power_curves.png
#
# Run:  Rscript code/14_Wobble_Analysis/11_power_sim_summarize.R
# ============================================================================

suppressPackageStartupMessages(library(tidyverse))

root <- rprojroot::find_root(rprojroot::is_git_root)
setwd(root)

simdir <- file.path("output_data", "14_Wobble_Analysis", "BAM_P0", "POWER_SIM")

design <- read_tsv(file.path(simdir, "sim_design.tsv"), show_col_types = FALSE)
res <- list.files(file.path(simdir, "results"), pattern = "^chunk_.*\\.tsv$", full.names = TRUE) %>%
  map_dfr(~ read_tsv(.x, col_names = c("sim_id", "p_wald"), show_col_types = FALSE))
stopifnot(!anyDuplicated(res$sim_id))
cat("sims with results:", nrow(res), "of", nrow(design), "\n")

thr <- read_tsv(file.path("output_data", "14_Wobble_Analysis", "BAM_P0", "SUMMARY", "scan_summary.tsv"),
                show_col_types = FALSE) %>%
  select(scan, thr_sig)
design_scan <- c(CLS = "CLS_BIN", MANUAL = "MANUAL_BIN", MANUAL_CLEAN = "MANUAL_BIN", QUANT = "VP0_INT")

d <- design %>%
  inner_join(res, by = "sim_id") %>%
  mutate(scan = design_scan[design]) %>%
  left_join(thr, by = "scan") %>%
  mutate(hit = -log10(p_wald) >= thr_sig)
stopifnot(!anyNA(d$thr_sig))

power <- d %>%
  group_by(design, scan, maf_bin, effect) %>%
  # GEMMA returns nan when the planted SNP separates cases completely (1 sim); dropped.
  summarise(n_rep = sum(!is.na(hit)), power = mean(hit, na.rm = TRUE), nominal_05 = mean(p_wald < 0.05, na.rm = TRUE),
            median_logp = median(-log10(p_wald), na.rm = TRUE), thr_sig = first(thr_sig), .groups = "drop")
write_tsv(power, file.path(simdir, "power_table.tsv"))

power80 <- power %>%
  group_by(design, maf_bin) %>%
  summarise(effect_80 = { x <- effect[power >= 0.8]; if (length(x)) min(x) else NA_real_ },
            max_power = max(power), .groups = "drop")
write_tsv(power80, file.path(simdir, "power_80.tsv"))

design_lab <- c(
  CLS          = "CLS: 16 v 42, classifier labels\n(x = per-allele odds ratio)",
  MANUAL       = "MANUAL: 25 v 41, 9 non-causal labels\n(x = per-allele odds ratio)",
  MANUAL_CLEAN = "MANUAL_CLEAN: 25 v 41, no label noise\n(x = per-allele odds ratio)",
  QUANT        = "QUANT: |vP0|, n = 66\n(x = variance explained by the SNP)"
)
plot_one <- function(dd, breaks, x_trans) {
  ggplot(dd, aes(effect, power, colour = maf_bin, shape = maf_bin)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
    annotate("text", x = max(dd$effect), y = 0.84, label = "80% power", hjust = 1, size = 3.2, colour = "grey40") +
    geom_line() + geom_point(size = 1.8) +
    facet_wrap(~ design, labeller = as_labeller(design_lab)) +
    scale_x_continuous(trans = x_trans, breaks = breaks) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = scales::percent,
                       expand = expansion(mult = c(0.02, 0.05))) +
    scale_colour_brewer(palette = "Dark2") +
    labs(x = NULL, y = "Power (genome-wide\npermutation threshold)", colour = "MAF bin", shape = "MAF bin") +
    theme_bw(base_size = 11) +
    theme(legend.text = element_text(size = 10))
}
p_bin <- plot_one(filter(power, design != "QUANT"), c(1, 2, 3, 5, 8, 12, 20, 50), "log")
p_q   <- plot_one(filter(power, design == "QUANT"), c(0, 0.1, 0.2, 0.3, 0.4), "identity")
p <- patchwork::wrap_plots(p_bin, p_q, ncol = 1, heights = c(1, 1), guides = "collect") +
  patchwork::plot_annotation(title = "Power of the BAM wobble scans to detect one locus")
ggsave(file.path(simdir, "power_curves.png"), p, width = 11, height = 8, dpi = 150)

cat("\nNull calibration (OR = 1 / PVE = 0):\n")
power %>% filter(effect %in% c(0, 1)) %>% select(design, maf_bin, power, nominal_05) %>% print()
cat("\nSmallest effect with >= 80% power:\n")
print(power80)
