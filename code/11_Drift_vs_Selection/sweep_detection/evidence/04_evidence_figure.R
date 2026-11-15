#!/usr/bin/env Rscript
## 04_evidence_figure.R
## Supplementary Fig. S7: the calibrated sweep evidence matrix (DESIGN.md sections 4-6).
##
## For each of the six tests (one per peak, over the lineage set in which its BP allele arose),
## the matched-null percentile of each primary statistic, the 95th and 99th percentile
## thresholds of the decision rule, and the resulting call with its empirical p.
##
## The unit of analysis is the one DESIGN.md section 2 defines. Peaks 2 and 3 are one
## test each over BP2+BP3 (the null takes the same maximum over the same two lineages), so the
## figure shows those rows and no per-lineage rows: evidence_matrix.tsv has none, and the
## per-lineage percentiles in peak_evidence_by_lineage.tsv are against the unmatched null (a),
## not the headline null, so they cannot sit on the same axis.
##
## All percentiles are oriented so that higher means more sweep-like: the pi ratio is a "low"
## statistic and 03_evidence_matrix.R already reverses its sense (pct = share of null regions
## with a ratio at least as high). The script checks that orientation, and recomputes the call,
## the composite score and the "2 of 3" rule from the percentiles, before plotting anything.
##
## The x axis is the upper-tail scale, -log10(1 - percentile/100), labelled in percentiles, so
## that the 95th-99th interval the rule turns on is readable. A percentile of 100 (no matched
## null region as extreme) is drawn at the null's resolution limit, 1 - 1/n.
##
## Usage (from the repo root; runs locally, needs data.table, ggplot2, patchwork, ragg):
##   Rscript code/11_Drift_vs_Selection/sweep_detection/evidence/04_evidence_figure.R [--in_dir DIR]
##
## Input:  evidence_matrix.tsv and calibration.tsv, from --in_dir if given, else from
##         output_data/11_Drift_vs_Selection/sweep_detection/evidence/ (where 03_evidence_matrix.R
##         writes them)
## Output: output_data/11_Drift_vs_Selection/sweep_detection/evidence/FigS7_evidence_matrix.pdf
##         output_data/11_Drift_vs_Selection/sweep_detection/evidence/FigS7_evidence_matrix.png (300 dpi)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(rprojroot)
})

## .git is a directory in a clone and a file in a linked worktree, and rprojroot's
## has_file() criterion matches only the latter. Accept either, so the script runs
## wherever the repository is checked out.
root <- find_root(has_file(".git") | has_dir(".git"))

args <- commandArgs(trailingOnly = TRUE)
arg_of <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

out_dir <- file.path(root, "output_data/11_Drift_vs_Selection/sweep_detection/evidence")
candidates <- c(arg_of("--in_dir", NA_character_), out_dir)
candidates <- candidates[!is.na(candidates)]
in_dir <- candidates[file.exists(file.path(candidates, "evidence_matrix.tsv"))][1]
if (is.na(in_dir)) stop("evidence_matrix.tsv not found in: ", paste(candidates, collapse = ", "))
cat("Reading evidence_matrix.tsv and calibration.tsv from\n  ", in_dir, "\n", sep = "")

em  <- fread(file.path(in_dir, "evidence_matrix.tsv"))
cal <- fread(file.path(in_dir, "calibration.tsv"))

PRIMARY     <- c("rsb", "clr", "pi_ratio")   # one per evidence class (H, S, D), DESIGN.md section 4
SUPPORT_PCT <- 95
STRONG_PCT  <- 99

## Peak genes and lineage labels as used in the paper
GENES <- c(`1` = "c6orf132", `2` = "itih3", `3` = "cebpg/cebpa",
           `4` = "cntnap5", `5` = "sptbn4b", `6` = "naa40")

## ------------------------------------------------------------------------------- checks
pm <- em[primary == TRUE]
stopifnot(setequal(unique(pm$statistic), PRIMARY),
          pm[, .N, by = peak]$N == length(PRIMARY),
          !anyNA(pm$pct_matched))

## Orientation: a low pi ratio is sweep-like, so the file must carry direction "low" for it and
## "high" for the other two. Within a peak the null is fixed, so this is what makes a high
## pct_matched mean "more sweep-like" for all three statistics.
stopifnot(all(pm[statistic == "pi_ratio", direction] == "low"),
          all(pm[statistic != "pi_ratio", direction] == "high"))

## Recompute the decision rule from the percentiles and refuse to plot a stale or edited matrix
rule <- pm[, .(n_ge95    = sum(pct_matched >= SUPPORT_PCT),
               composite_re = sort(pct_matched, decreasing = TRUE)[2],
               any_ge99  = any(pct_matched >= STRONG_PCT)),
           by = .(peak, lineages, composite, supported, strong, empirical_p, tp_label,
                  n_matched_null)]
rule[, `:=`(supported_re = n_ge95 >= 2, strong_re = n_ge95 >= 2 & any_ge99)]
bad <- rule[abs(composite_re - composite) > 1e-9 | supported_re != supported | strong_re != strong]
if (nrow(bad)) {
  print(bad)
  stop("the call or composite in evidence_matrix.tsv does not follow from its percentiles")
}
rule <- merge(rule, cal[, .(peak, fpr_supported, fpr_strong)], by = "peak", all.x = TRUE)
setorder(rule, peak)

## ------------------------------------------------------------------------ plotting data
tail_x <- function(pct) -log10(1 - pct / 100)   # upper-tail scale; 95 -> 1.30, 99 -> 2.00

STAT_LEVELS <- c("rsb", "clr", "pi_ratio")
STAT_LABELS <- c(rsb = "Rsb (haplotype)", clr = "SweepFinder2 CLR (site frequency)",
                 pi_ratio = "π(BP) / π(TP) (diversity)")
STAT_OFFSET <- c(rsb = 0.27, clr = 0, pi_ratio = -0.27)   # fixed sub-row order within a peak

pm[, stat := factor(statistic, levels = STAT_LEVELS)]
pm[, y := -(peak - 1) + STAT_OFFSET[statistic]]
## a percentile of 100 has no finite tail value; draw it at the resolution of its null (1/n)
pm[, at_limit := pct_matched >= 100 - 1e-9]
pm[, x := tail_x(pmin(pct_matched, 100 * (1 - 1 / n_matched_null)))]
## direct labels only on the points the rule counts (>= 95th); "100*" flags a point drawn at
## the resolution limit rather than at its own value (the figure legend explains the asterisk)
pm[, lab := ifelse(pct_matched >= SUPPORT_PCT,
                   ifelse(at_limit, "100*", formatC(pct_matched, format = "f", digits = 1)), NA)]

rule[, y := -(peak - 1)]
rule[, call := fifelse(strong, "Strong", fifelse(supported, "Supported", "Not supported"))]
rule[, p_lab := formatC(signif(empirical_p, 2), format = "fg", digits = 2, flag = "#")]  # keeps 0.060
rule[, lineage_lab := gsub("+", " + ", lineages, fixed = TRUE)]

cat("\nPlotted values (matched-null percentile; call; empirical p):\n")
print(dcast(pm, peak + lineages ~ statistic, value.var = "pct_matched")[
  rule[, .(peak, call, empirical_p, n_matched_null, fpr_supported, fpr_strong)], on = "peak"])

## ---------------------------------------------------------------------------- tokens
INK      <- "#0b0b0b"
INK_2    <- "#45443f"
GRID     <- "#e1e0d9"
BAND_95  <- "#f3f2ee"
BAND_99  <- "#dcdbd3"   # two greys far enough apart to survive greyscale print
RULE_COL <- "#898781"
SERIES   <- c(rsb = "#2a78d6", clr = "#eb6834", pi_ratio = "#1baf7a")  # validated, 3 slots all-pairs
SHAPES   <- c(rsb = 21, clr = 24, pi_ratio = 22)                        # redundant with colour

FAMILY <- "Arial"
BASE   <- 7.5                       # pt, at the 180-mm print width
pt     <- function(x) x / .pt       # ggplot text size in mm from pt

y_top    <- 0.80                    # header row (column titles, threshold labels)
y_lim    <- c(-5.55, 1.0)
x_lim    <- c(tail_x(50) - 0.12, tail_x(99.97))   # room for the "100" label
x_breaks <- tail_x(c(50, 75, 90, 95, 99, 99.9))
x_labels <- c("50", "75", "90", "95", "99", "99.9")
seps     <- -(1:5) + 0.5            # hairlines between peaks

theme_panel <- theme_void(base_family = FAMILY, base_size = BASE) +
  theme(plot.margin = margin(0, 0, 0, 0))

## ------------------------------------------------------------------- left: peak labels
left <- ggplot(rule) +
  geom_text(aes(x = 0, y = y + 0.15,
                label = sprintf("bold('Peak %d  ')*italic('%s')", peak, GENES[as.character(peak)])),
            parse = TRUE, hjust = 0, family = FAMILY, size = pt(BASE), colour = INK) +
  geom_text(aes(x = 0, y = y - 0.2, label = lineage_lab),
            hjust = 0, family = FAMILY, size = pt(BASE), colour = INK_2) +
  geom_hline(yintercept = seps, colour = GRID, linewidth = 0.25) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
  theme_panel

## --------------------------------------------------------- middle: percentiles vs thresholds
mid <- ggplot(pm) +
  annotate("rect", xmin = tail_x(SUPPORT_PCT), xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = BAND_95) +
  annotate("rect", xmin = tail_x(STRONG_PCT), xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = BAND_99) +
  geom_vline(xintercept = tail_x(c(50, 75, 90)), colour = GRID, linewidth = 0.25) +
  geom_vline(xintercept = tail_x(c(SUPPORT_PCT, STRONG_PCT)), colour = RULE_COL, linewidth = 0.35) +
  geom_hline(yintercept = seps, colour = GRID, linewidth = 0.25) +
  ## plain "95th"/"99th": cairo_pdf drew a "≥" glyph from a fallback font
  annotate("text", x = tail_x(SUPPORT_PCT) + 0.04, y = y_top, label = "95th",
           hjust = 0, family = FAMILY, size = pt(BASE), colour = INK_2) +
  annotate("text", x = tail_x(STRONG_PCT) + 0.04, y = y_top, label = "99th",
           hjust = 0, family = FAMILY, size = pt(BASE), colour = INK_2) +
  geom_point(aes(x = x, y = y, fill = stat, shape = stat),
             colour = "white", stroke = 0.4, size = 2.4) +
  geom_text(aes(x = x + 0.07, y = y, label = lab), na.rm = TRUE,
            hjust = 0, family = FAMILY, size = pt(BASE), colour = INK_2) +
  scale_fill_manual(values = SERIES, labels = STAT_LABELS, name = NULL) +
  scale_shape_manual(values = SHAPES, labels = STAT_LABELS, name = NULL) +
  scale_x_continuous(limits = x_lim, breaks = x_breaks, labels = x_labels, expand = c(0, 0),
                     name = paste0("Percentile of matched null (higher = more sweep-like; π ratio reversed)\n",
                                   "log scale of 100 – percentile")) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
  guides(fill = guide_legend(override.aes = list(size = 2.6, colour = "white"))) +
  theme_minimal(base_family = FAMILY, base_size = BASE) +
  theme(panel.grid      = element_blank(),
        axis.text.y     = element_blank(),
        axis.title.y    = element_blank(),
        axis.text.x     = element_text(colour = INK_2, size = BASE),
        axis.title.x    = element_text(colour = INK, size = BASE, margin = margin(t = 3),
                                       lineheight = 1.1),
        axis.line.x     = element_line(colour = "#c3c2b7", linewidth = 0.3),
        axis.ticks.x    = element_line(colour = "#c3c2b7", linewidth = 0.3),
        axis.ticks.length.x = unit(1.2, "mm"),
        legend.position = "top",
        legend.justification = "left",
        legend.text     = element_text(colour = INK, size = BASE),
        legend.key.size = unit(3.2, "mm"),
        legend.spacing.x = unit(1, "mm"),
        legend.margin   = margin(0, 0, 1, 0),
        legend.box.spacing = unit(1, "mm"),
        plot.margin     = margin(0, 2, 0, 2))

## ------------------------------------------------------------- right: call and empirical p
right <- ggplot(rule) +
  ## both headers go through plotmath, and the phantom p gives "Call" the same descender, so
  ## the two share a baseline. Spaces stay inside the quoted strings, here and in the peak
  ## labels: a plotmath "~" made cairo_pdf embed a fallback font for the space glyph.
  annotate("text", x = 0.04, y = y_top, label = "bold('Call')*phantom(bolditalic(p))", parse = TRUE,
           hjust = 0, family = FAMILY, size = pt(BASE), colour = INK) +
  annotate("text", x = 0.56, y = y_top, label = "bold('Empirical ')*bolditalic(p)", parse = TRUE,
           hjust = 0, family = FAMILY, size = pt(BASE), colour = INK) +
  geom_text(aes(x = 0.04, y = y, label = call, fontface = ifelse(supported, "bold", "plain"),
                colour = supported),
            hjust = 0, family = FAMILY, size = pt(BASE)) +
  geom_text(aes(x = 0.56, y = y, label = p_lab, fontface = ifelse(supported, "bold", "plain"),
                colour = supported),
            hjust = 0, family = FAMILY, size = pt(BASE)) +
  geom_hline(yintercept = seps, colour = GRID, linewidth = 0.25) +
  scale_colour_manual(values = c(`TRUE` = INK, `FALSE` = INK_2), guide = "none") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
  theme_panel

## ------------------------------------------------------------------------------ layout
W_MM <- 180; H_MM <- 88
fig <- left + mid + right +
  plot_layout(widths = c(30, 108, 42), guides = "collect") &
  theme(legend.position = "top", legend.justification = "left",
        plot.background = element_rect(fill = "white", colour = NA))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pdf_out <- file.path(out_dir, "FigS7_evidence_matrix.pdf")
png_out <- file.path(out_dir, "FigS7_evidence_matrix.png")
ggsave(pdf_out, fig, width = W_MM, height = H_MM, units = "mm", device = cairo_pdf)
ggsave(png_out, fig, width = W_MM, height = H_MM, units = "mm", dpi = 300, device = ragg::agg_png,
       bg = "white")
cat("\nWrote (", W_MM, " x ", H_MM, " mm)\n  ", pdf_out, "\n  ", png_out, "\n", sep = "")
