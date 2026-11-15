## 05_summarise.R
## Peak 6 CURV3 summary and verdict (ESM, Targeted test at Peak 6).
##   - Dtrios: one trio per test, signed so + = sharing between the named lineage and CURV3 (orient_trios, as in
##     DsuiteAnalysis.Rmd). D_genome = DtriosCombine over chr1-chr25 for the same trio; Z_vs_genome =
##     (D_peak - D_genome) / SE_peak with SE_peak = |D/Z| of the peak run, as Data S7. q = BH over the six tests.
##   - f_dM: Dinvestigate -w 50,25 windows. Peak median of the windows overlapping the extended Peak 6 interval, compared
##     with all background windows and with medians of equal-size blocks of background windows (windows overlapping any
##     of the six peaks excluded). The code is the main analysis's: summarize_rep_side_globalbg() from
##     dinvestigate_functions.R and build_peak_summary() read from DsuiteInvestigateAnalysis.Rmd at run time. Trios are
##     P1 = tested lineage, so f_dM is multiplied by -1 (sign_flip = -1, the main analysis's representative side).
##   - Rule (as for the main tests, plus the f_dM part): signal = q < 0.05 & D > 0 & Z_vs_genome > 2 & f_dM peak median
##     above the 95th block percentile. D_signal is the D part alone (the definition of the `signal` column in Data S7).
##   - Verdict branches. "The TP controls show the same excess" is read as: TP1-CURV3 or TP2-CURV3 meets the same full
##     rule; the verdict text also says if a TP control meets only the D part.
## Also checks the summary code against the main analysis (peak_fdm_summary.tsv, BP2-BP_ANC Peak 6) and the chr24
## Dinvestigate rerun against output_data/10_Introgression/BP2_TP_OUT_BP_ANC_localFstats_chr24_50_25.txt.

suppressPackageStartupMessages(library(tidyverse))

root <- system("git rev-parse --show-toplevel", intern = TRUE)
intro_dir <- file.path(root, "output_data", "10_Introgression")
out_dir   <- file.path(intro_dir, "peak6_curv3")

source(file.path(root, "code/10_Introgression/dinvestigate_functions.R"))
rmd <- readLines(file.path(root, "code/10_Introgression/DsuiteInvestigateAnalysis.Rmd"))
i0 <- grep("^build_peak_summary <- function", rmd)
i1 <- i0 + which(rmd[(i0 + 1):length(rmd)] == "}")[1]
eval(parse(text = rmd[i0:i1]))

## ---------------- Dtrios ----------------
read_dsuite <- function(f, peak_label = NA_character_) {
  df <- read_tsv(f, show_col_types = FALSE, progress = FALSE)
  names(df)[1:10] <- c("P1","P2","P3","D","Z","p","f4ratio","BBAA","ABBA","BABA")
  df %>% mutate(Peak = peak_label)
}

hyp_defs <- tribble(
  ~Hypothesis,             ~focal, ~sister,  ~P3,     ~type,
  "BP2–CURV3",             "BP2",  "TP_OUT", "CURV3", "primary",
  "BP2–CURV3 (P1 = TP2)",  "BP2",  "TP2",    "CURV3", "sister form",
  "TP1–CURV3",             "TP1",  "TP_OUT", "CURV3", "TP control",
  "TP2–CURV3",             "TP2",  "TP_OUT", "CURV3", "TP control",
  "BP1–CURV3",             "BP1",  "TP_OUT", "CURV3", "BP control",
  "BP3–CURV3",             "BP3",  "TP_OUT", "CURV3", "BP control"
)

orient_trios <- function(df) {
  hyp_defs %>%
    pmap_dfr(function(Hypothesis, focal, sister, P3, type) {
      hit <- df[df$P3 == P3 & ((df$P1 == focal & df$P2 == sister) | (df$P1 == sister & df$P2 == focal)), ]
      stopifnot(nrow(hit) == n_distinct(df$Peak))
      hit %>% mutate(Hypothesis = Hypothesis, type = type, focal = focal, sister = sister,
                     sign = ifelse(P2 == focal, 1, -1))
    }) %>%
    mutate(D_o = sign * D, Z_o = sign * Z, f4_o = sign * f4ratio,
           SE = ifelse(Z != 0, abs(D / Z), NA_real_),
           nSites = ABBA + BABA)
}

peak_raw <- read_dsuite(file.path(out_dir, "peak_curv3", "peak_6_tree.txt"), "Peak 6")
gw_raw   <- read_dsuite(file.path(out_dir, "genome_curv3", "whole_genome_combined_tree.txt"), "genome")
# The tree must put CURV3 as P3 in every trio that contains it
stopifnot(all(with(peak_raw, P3[P1 == "CURV3" | P2 == "CURV3" | P3 == "CURV3"]) == "CURV3"),
          all(with(gw_raw,   P3[P1 == "CURV3" | P2 == "CURV3" | P3 == "CURV3"]) == "CURV3"))

gw <- gw_raw %>% orient_trios() %>% select(Hypothesis, D_gw = D_o, Z_gw = Z_o)
dt <- peak_raw %>% orient_trios() %>%
  left_join(gw, by = "Hypothesis") %>%
  mutate(Z_diff = (D_o - D_gw) / SE,
         q = p.adjust(p, method = "BH"))
stopifnot(nrow(dt) == 6)

## ---------------- f_dM (Dinvestigate) ----------------
read_dinvestigate <- function(f, peak_label = NA_character_) {
  df <- read_tsv(f, show_col_types = FALSE, progress = FALSE)
  names(df) <- c("chr", "windowStart", "windowEnd", "D", "f_d", "f_dM", "d_f")
  df %>% mutate(Comparision = peak_label)
}
load_dinv <- function(dir, comps) {
  files <- list.files(dir, pattern = "_localFstats_chr[0-9]+_50_25\\.txt(\\.gz)?$", full.names = TRUE)
  lab <- str_remove(basename(files), "_localFstats_.*$")
  files <- files[lab %in% comps]; lab <- lab[lab %in% comps]
  stopifnot(all(table(lab)[comps] == 25))
  map2_dfr(files, lab, read_dinvestigate)
}

peaks.bed <- read.table(file.path(root, "input_data/10_Introgression/peaks.bed"))
names(peaks.bed) <- c("chrom", "start", "end", "peak", "idx_p", "idx_bp")
peaks.bed$chrom <- paste0("chr", peaks.bed$chrom)
peaks.bed <- peaks.bed[order(peaks.bed$peak), ]

fdm_summary <- function(dinv_raw, rep_map) {
  all_summ <- summarize_rep_side_globalbg(dinv_raw, peaks.bed, rep_map, bg_stat = median)
  build_peak_summary(dinv_raw, peaks.bed, rep_map, all_summ)
}
make_rep_map_fixed <- function(defs) {
  defs %>% transmute(Hypothesis, Comparision = paste(focal, sister, P3, sep = "_"), P1 = focal, Spacer = sister, P3,
                     trio_label = paste0("(", P1, ", ", Spacer, ")(", P3, ")"), sign_flip = -1)
}

# Check 1: the same code on the main analysis's BP2-TP_OUT-BP_ANC windows must give peak_fdm_summary.tsv (BP2↔BP_ANC)
pub_map  <- tibble(Hypothesis = "BP2↔BP_ANC", focal = "BP2", sister = "TP_OUT", P3 = "BP_ANC") %>% make_rep_map_fixed()
pub_fdm  <- fdm_summary(load_dinv(intro_dir, "BP2_TP_OUT_BP_ANC"), pub_map)
pub_ref  <- read_tsv(file.path(intro_dir, "peak_fdm_summary.tsv"), show_col_types = FALSE) %>% filter(Hypothesis == "BP2↔BP_ANC")
check_fdm <- inner_join(pub_ref, pub_fdm, by = c("Hypothesis", "Peaks"), suffix = c("_published", "_rerun"))
stopifnot(nrow(check_fdm) == 6,
          isTRUE(all.equal(check_fdm$df_median_published, check_fdm$df_median_rerun)),
          isTRUE(all.equal(check_fdm$percentile_block_published, check_fdm$percentile_block_rerun)),
          all(check_fdm$nwin_published == check_fdm$nwin_rerun))

# Check 2: chr24 Dinvestigate rerun of the main analysis's trio vs its published window file
f_old <- file.path(intro_dir, "BP2_TP_OUT_BP_ANC_localFstats_chr24_50_25.txt")
f_new <- list.files(file.path(out_dir, "dinvestigate_published_check"), pattern = "^BP2_TP_OUT_BP_ANC_localFstats_chr24_50_25\\.txt(\\.gz)?$",
                    full.names = TRUE)
w_old <- read_dinvestigate(f_old); w_new <- read_dinvestigate(f_new)
check_win <- tibble(
  file = basename(f_old), windows_published = nrow(w_old), windows_rerun = nrow(w_new),
  identical_coordinates = nrow(w_old) == nrow(w_new) && all(w_old$windowStart == w_new$windowStart & w_old$windowEnd == w_new$windowEnd),
  max_abs_diff_fdM = if (nrow(w_old) == nrow(w_new)) max(abs(w_old$f_dM - w_new$f_dM), na.rm = TRUE) else NA_real_,
  max_abs_diff_D   = if (nrow(w_old) == nrow(w_new)) max(abs(w_old$D - w_new$D), na.rm = TRUE) else NA_real_)
write_tsv(check_fdm %>% select(Hypothesis, Peaks, nwin_published, nwin_rerun, df_median_published, df_median_rerun,
                               percentile_block_published, percentile_block_rerun),
          file.path(out_dir, "check_fdm_code_vs_published_peak_fdm_summary.tsv"))
write_tsv(check_win, file.path(out_dir, "check_dinvestigate_chr24_vs_published.tsv"))
print(as.data.frame(check_win))

# The six CURV3 trios
rep_map <- make_rep_map_fixed(hyp_defs)
dinv_raw <- load_dinv(file.path(out_dir, "dinvestigate_curv3"), rep_map$Comparision)
fdm6 <- fdm_summary(dinv_raw, rep_map) %>% filter(Peaks == 6)
stopifnot(nrow(fdm6) == 6)

## ---------------- Table and rule ----------------
summ <- dt %>%
  left_join(fdm6 %>% select(Hypothesis, df_median, nwin, bg_median, percentile, percentile_block), by = "Hypothesis") %>%
  mutate(D_signal = q < 0.05 & D_o > 0 & Z_diff > 2,
         signal = D_signal & percentile_block > 95,
         Hypothesis = factor(Hypothesis, levels = hyp_defs$Hypothesis)) %>%
  arrange(Hypothesis)

tab <- summ %>%
  transmute(Test = Hypothesis, type, P1, P2, P3,
            D = round(D_o, 3), Z = round(Z_o, 2), p = signif(p, 2), q = signif(q, 2),
            f4_ratio = round(f4_o, 3), D_genome = round(D_gw, 3), Z_vs_genome = round(Z_diff, 2),
            informative_sites = round(nSites),
            fdM_peak_median = round(df_median, 3), fdM_nwin = nwin, fdM_bg_median = round(bg_median, 3),
            fdM_window_percentile = round(percentile, 1), fdM_block_percentile = round(percentile_block, 1),
            D_signal, signal)
write_tsv(tab, file.path(out_dir, "peak6_curv3_dtrios_summary.tsv"))
print(as.data.frame(tab))

## ---------------- Verdict ----------------
row_of <- function(h) summ %>% filter(Hypothesis == h)
prim <- row_of("BP2–CURV3"); sis <- row_of("BP2–CURV3 (P1 = TP2)")
tp   <- summ %>% filter(type == "TP control")
tp_full <- any(tp$signal); tp_Donly <- any(tp$D_signal & !tp$signal)

fmt <- function(r) sprintf("D = %.3f, Z = %.2f, q = %.2g, %.2f SE vs genome-wide D = %.3f, f_dM median %.3f at the %.1fth block percentile",
                           r$D_o, r$Z_o, r$q, r$Z_diff, r$D_gw, r$df_median, r$percentile_block)

if (prim$signal && !tp_full) {
  branch <- "Branch 1 (introgression signal): BP2-CURV3 meets the rule and neither TP control does."
} else if (prim$signal && tp_full) {
  branch <- "Branch 2 (shared ancestry, not BP-specific): BP2-CURV3 meets the rule, but a TP control shows the same excess."
} else {
  branch <- "Branch 3 (cannot be told apart from shared ancestral variation): the primary test BP2-CURV3 does not meet the rule."
}
crit_failed <- function(r) {
  f <- c("q < 0.05"[!(r$q < 0.05)], "D > 0"[!(r$D_o > 0)], "> 2 SE above genome-wide D"[!(r$Z_diff > 2)],
         "f_dM above the 95th block percentile"[!(r$percentile_block > 95)])
  if (length(f)) paste0(" (fails: ", paste(f, collapse = ", "), ")") else ""
}
s2 <- sprintf("Primary BP2-CURV3 (P1 = TP_OUT) at Peak 6: %s%s; TP controls: TP1-CURV3 rule %s (D part %s), TP2-CURV3 rule %s (D part %s)%s.",
              fmt(prim), crit_failed(prim),
              ifelse(tp$signal[tp$Hypothesis == "TP1–CURV3"], "met", "not met"),
              ifelse(tp$D_signal[tp$Hypothesis == "TP1–CURV3"], "met", "not met"),
              ifelse(tp$signal[tp$Hypothesis == "TP2–CURV3"], "met", "not met"),
              ifelse(tp$D_signal[tp$Hypothesis == "TP2–CURV3"], "met", "not met"),
              ifelse(tp_Donly, "; a TP control meets the D part of the rule but not the f_dM part", ""))
# The TP_OUT trio is the primary test and the sister form is not used to set the branch. Say so whenever the
# sister form disagrees with the primary.
s3 <- sprintf("Sister form (P1 = TP2): %s%s, rule %s%s",
              fmt(sis), crit_failed(sis), ifelse(sis$signal, "met", "not met"),
              ifelse(sis$signal != prim$signal,
                     "; the TP_OUT trio is the primary test, so the branch above rests on it and the sister-form result is not used to change it.",
                     "."))
verdict <- paste(branch, s2, s3)
writeLines(c(verdict, "",
             "Rule: signal = BH q < 0.05 across the six tests, D > 0, D more than 2 SE above the genome-wide D for the trio, and f_dM peak median above the 95th block percentile. Branches: signal and no TP-control signal = introgression signal; TP controls show the same excess = shared ancestry, not BP-specific; otherwise the sharing cannot be told apart from shared ancestral variation."),
           file.path(out_dir, "rule_verdict.txt"))
cat("\n", verdict, "\n")
