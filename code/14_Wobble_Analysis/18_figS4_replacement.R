#!/usr/bin/env Rscript
# ============================================================================
# 18_figS4_replacement.R
#
# Fig. S4 (fine-scale haplotypes at the six GWAS peaks).
# The SNP set (the TOP250 VCF; see README) and the population tree are those of
# code/08_Peak_Analysis/phylogenoplot_multi.R (called from
# peak_analysis_new.Rmd), but each row is one PHASED haplotype and each cell is
# coloured by BP vs TP allele instead of the unphased REF/HET/ALT genotype:
#   - Peaks 1, 4, 5: BP allele = Bambomo BP consensus allele.
#   - Peaks 2, 3: BP allele = Bikagala/Bavavela/Biroundou BP consensus allele.
#   - Peak 6 (BP2-only locus): BP allele = Bikagala/Bavavela BP consensus allele.
#   - Only diagnostic SNPs (BP and Apassa/Bengue/Douvalou/Douengi TP consensus
#     differ) are drawn.
# Each fish has two rows; within a fish, the more BP-like copy is drawn on top
# (phase is not comparable across chromosomes). Fish keep the same order in
# every panel: population (tree order), then phenotype (BP, mixed, TP), then ID.
# Orange diamonds right of a panel mark copies, outside the peak's home
# populations, that are the SAME haplotype as a home BP-class copy by the
# identity test of 19_haplotype_identity_test.R (Fig. S4 SNPs): identical, or
# 1-2 mismatches at alleles that home copies also carry. Similar but
# distinguishable copies are not marked. The green bar marks the SNPs inside
# the current peak interval. Tip labels give the number of fish.
#
# Population assignment follows the sample metadata (POP), not the ID prefix,
# so APA_223 (collected at Bambomo Creek) is drawn with Bambomo.
#
# Inputs:
#   output_data/14_Wobble_Analysis/PKINGS_WOB_MULTILOCUS/peak_hap_zoom/figs4_haplotypes.rds
#     (written by 17_peak_hap_match_audit.R)
#   .../peak_hap_zoom/identity_test_figS4snps.tsv (written by 19_haplotype_identity_test.R)
#   input_data/04_Phylogenomics/pkings_iqtree_genes_astral_gcf_scf_tree.cf.tree.ultra.time.nex
#
# Outputs (PKINGS_WOB_MULTILOCUS/peak_hap_zoom/):
#   figS4_replacement.png / .pdf
#
# Run:  Rscript code/14_Wobble_Analysis/18_figS4_replacement.R
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ape)
  library(ggtree)
  library(patchwork)
})

root <- rprojroot::find_root(rprojroot::has_file(".git") | rprojroot::is_git_root)
setwd(root)
outdir <- file.path(root, "output_data", "14_Wobble_Analysis", "PKINGS_WOB_MULTILOCUS", "peak_hap_zoom")

dat <- readRDS(file.path(outdir, "figs4_haplotypes.rds"))
figs4 <- dat$figs4; figs4_long <- dat$figs4_long; bed <- dat$bed; key <- dat$key
# Same haplotype = identical to a home copy, or 1-2 mismatches with no allele
# that is absent from every home copy (19_haplotype_identity_test.R).
same_ids <- read_tsv(file.path(outdir, "identity_test_figS4snps.tsv"), show_col_types = FALSE) %>%
  filter(nn_mismatch == 0 | (nn_mismatch <= 2 & novel_alleles == 0)) %>%
  transmute(peak, id)

GAP <- 4        # blank rows between populations
HEAD <- 8       # rows of headroom above the top population (peak-interval bar)
home_pops <- list(BAM = "BAM", other = c("BIK", "BAVA", "BIR"))

# ---------------------------------------------------------------------------
# 1. Population tree (as in peak_analysis_new.Rmd), tips renamed to POP codes
# ---------------------------------------------------------------------------
tree <- read.tree(file.path(root, "input_data", "04_Phylogenomics",
                            "pkings_iqtree_genes_astral_gcf_scf_tree.cf.tree.ultra.time.nex"))
tree <- drop.tip(tree, c("BIK_6927", "COB_4006", "PHOP_5511"))
tree$tip.label <- sub("_.*$", "", tree$tip.label)
stopifnot(!anyDuplicated(tree$tip.label), all(unique(figs4$POP) %in% tree$tip.label))

# ---------------------------------------------------------------------------
# 2. Row layout: populations in tree order (top to bottom), fish within
# ---------------------------------------------------------------------------
# Not ladderized: populations run top-down in newick order (BAM, APA, BEN, ...),
# the order of phylogenoplot_multi.R. ggtree puts the first tip at the bottom.
tip_y <- ggtree(tree, ladderize = FALSE)$data %>% filter(isTip) %>% select(label, y)
pop_order <- tip_y %>% arrange(y) %>% pull(label)

pheno_lab <- function(Ph) case_when(Ph == "BP" ~ "BP", Ph == "Wobble" ~ "Mixed", TRUE ~ "TP")
fish_rows <- figs4 %>%
  distinct(PN, POP, Phenotype) %>%
  mutate(Pheno = factor(pheno_lab(Phenotype), levels = c("BP", "Mixed", "TP")),
         POP = factor(POP, levels = pop_order)) %>%
  arrange(POP, Pheno, PN)

# Top-down y: each fish takes two rows; GAP rows between populations.
pop_sizes <- fish_rows %>% count(POP) %>% mutate(rows = 2 * n)
pop_sizes$top <- NA_real_
y_cursor <- 0
for (i in seq_len(nrow(pop_sizes))) {
  pop_sizes$top[i] <- y_cursor
  y_cursor <- y_cursor + pop_sizes$rows[i] + GAP
}
total_rows <- y_cursor - GAP
fish_rows <- fish_rows %>%
  group_by(POP) %>% mutate(i = row_number()) %>% ungroup() %>%
  left_join(pop_sizes %>% select(POP, top), by = "POP") %>%
  mutate(y1 = -(top + 2 * (i - 1) + 1), y2 = y1 - 1)   # negative: top of plot = row 1
pop_center <- pop_sizes %>% mutate(center = -(top + (rows + 1) / 2))
y_lim <- c(-(total_rows + 3), HEAD)

# ---------------------------------------------------------------------------
# 3. Tree panel: tips at their population's block centre
# ---------------------------------------------------------------------------
n_fish <- fish_rows %>% count(POP)
p_tree <- ggtree(tree, linewidth = 0.55, ladderize = FALSE)
newy <- numeric(nrow(p_tree$data))
tips <- p_tree$data$isTip
newy[tips] <- pop_center$center[match(p_tree$data$label[tips], as.character(pop_center$POP))]
# Internal nodes: mean of their children's y, deepest first. Children come
# from ggtree's own parent column (its node numbering).
td <- p_tree$data
for (nd in td$node[!tips][order(-td$x[!tips])]) {
  kids <- td$node[td$parent == nd & td$node != nd]
  newy[td$node == nd] <- mean(newy[match(kids, td$node)])
}
p_tree$data$y <- newy
p_tree$data$tiplab <- ifelse(p_tree$data$isTip,
  sprintf("%s (%d)", p_tree$data$label, n_fish$n[match(p_tree$data$label, as.character(n_fish$POP))]), NA)
x_max <- max(p_tree$data$x)
p_tree <- p_tree +
  geom_tiplab(aes(label = tiplab), size = 2.7, offset = x_max * 0.03) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.42))) +
  theme_tree() +
  theme(plot.margin = margin(2, 0, 2, 2))

# ---------------------------------------------------------------------------
# 4. Phenotype strip (Okabe-Ito colours)
# ---------------------------------------------------------------------------
pheno_cols <- c(BP = "#0072B2", Mixed = "#CC79A7", TP = "#D55E00")
strip_d <- bind_rows(fish_rows %>% transmute(Pheno, y = y1), fish_rows %>% transmute(Pheno, y = y2))
p_strip <- ggplot(strip_d, aes(x = 1, y = y, fill = Pheno)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = pheno_cols, name = "Phenotype") +
  scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(plot.margin = margin(2, 4, 2, 0))

# ---------------------------------------------------------------------------
# 5. Peak panels
# ---------------------------------------------------------------------------
col_bp   <- "#262626"   # BP allele
col_tp   <- "#E0E0E0"   # TP allele
col_peak <- "#009E73"   # current peak interval
col_flag <- "#E69F00"   # same haplotype as a home BP copy, outside the home populations
key_peak <- "SNPs in the current peak interval"
key_flag <- "Copy outside the home populations that is the\nsame haplotype as a home BP copy (identity test)"

peak_panel <- function(pk, colour_guide) {
  hs <- figs4 %>% filter(peak == pk) %>%
    group_by(PN) %>% arrange(desc(bp_match_frac), hap, .by_group = TRUE) %>%
    mutate(copy = row_number()) %>% ungroup() %>%
    left_join(fish_rows %>% select(PN, y1, y2), by = "PN") %>%
    mutate(y = if_else(copy == 1, y1, y2),
           home_pop = (origin == "BAM" & POP %in% home_pops$BAM) |
                      (origin == "other" & POP %in% if (pk == 6) c("BIK", "BAVA") else home_pops$other))
  d <- figs4_long %>% filter(peak == pk) %>% inner_join(hs %>% select(id, y), by = "id")
  n_diag <- max(d$snp)
  cur <- bed[bed$peak == pk, ]
  cur_snps <- range(d$snp[d$pos >= cur$start & d$pos <= cur$end])
  flags <- hs %>% filter(id %in% same_ids$id[same_ids$peak == pk], !home_pop) %>%
    mutate(x = n_diag + n_diag * 0.06, key = key_flag)
  bar <- tibble(x = cur_snps[1] - 0.5, xend = cur_snps[2] + 0.5, y = 3.5, key = key_peak)
  origin_lab <- if (key$origin[key$peak == pk] == "BAM") "Bambomo" else if (pk == 6) "BIK/BAVA" else "BIK/BAVA/BIR"
  flag_lab <- if (nrow(flags) == 0) "same haplotype: none" else
    sprintf("same haplotype: %d cop%s, %d fish", nrow(flags), if (nrow(flags) == 1) "y" else "ies",
            n_distinct(flags$PN))
  ggplot(d, aes(snp, y)) +
    geom_tile(aes(fill = bp), width = 1, height = 1) +
    geom_segment(data = bar, aes(x = x, xend = xend, y = y, yend = y, colour = key), linewidth = 2.2) +
    geom_point(data = flags, aes(x = x, y = y, colour = key), shape = 18, size = 1.9) +
    scale_fill_manual(values = c(`TRUE` = col_bp, `FALSE` = col_tp),
                      labels = c(`TRUE` = "BP allele", `FALSE` = "TP allele"), name = "Allele") +
    scale_colour_manual(values = setNames(c(col_peak, col_flag), c(key_peak, key_flag)), name = NULL,
                        limits = c(key_peak, key_flag), drop = FALSE,
                        guide = guide_legend(override.aes = list(shape = c(NA, 18), linewidth = c(2.2, 0),
                                                                 size = c(NA, 3)))) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    scale_x_continuous(limits = c(0.5, n_diag + n_diag * 0.11), expand = c(0, 0),
                       breaks = function(l) { b <- pretty(c(1, n_diag), 3); b[b >= 1 & b <= n_diag] }) +
    labs(title = sprintf("Peak %d", pk),
         subtitle = sprintf("%s, %d SNPs\nBP = %s\n%s", key$chr[key$peak == pk], n_diag, origin_lab, flag_lab),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 8.5) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = 7),
          plot.title = element_text(size = 10.5, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 7.5, hjust = 0.5, lineheight = 0.95),
          plot.margin = margin(2, 4, 2, 4)) +
    (if (!colour_guide) guides(colour = "none") else NULL)
}

# The bar/flag legend is drawn once (patchwork does not merge it across panels).
panels <- lapply(1:6, function(pk) peak_panel(pk, colour_guide = pk == 1))

p <- wrap_plots(c(list(p_tree, p_strip), panels), nrow = 1,
                widths = c(1.5, 0.08, rep(1, 6))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.justification = "center",
        legend.text = element_text(size = 8), legend.title = element_text(size = 8.5),
        legend.margin = margin(6, 6, 0, 6))

ggsave(file.path(outdir, "figS4_replacement.png"), p, width = 11, height = 11, dpi = 400, bg = "white")
ggsave(file.path(outdir, "figS4_replacement.pdf"), p, width = 11, height = 11, bg = "white")
cat("Wrote", file.path(outdir, "figS4_replacement.png"), "\n")
