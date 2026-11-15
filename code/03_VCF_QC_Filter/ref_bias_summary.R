## ref_bias_summary.R
## Reduce the per-chromosome reference-bias tables (ref_bias_stats.sb) and compare BP vs TP fish.
## The reference is a BP (Bambomo) male; do TP reads map worse, especially in variable non-coding regions?
## Usage: Rscript code/03_VCF_QC_Filter/ref_bias_summary.R [ref_bias_dir]

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
root <- system("git rev-parse --show-toplevel", intern = TRUE)
args <- commandArgs(TRUE)
dir <- if (length(args)) args[1] else file.path(root, "output_data/03_QC/ref_bias")

raw <- bind_rows(lapply(list.files(file.path(dir, "per_chrom"), full.names = TRUE), read.delim))
pooled_peaks <- raw %>% filter(grepl("^peak", label)) %>% mutate(label = "peaks_all")
d <- bind_rows(raw, pooled_peaks) %>%
  group_by(label, sample) %>%
  summarise(across(c(n_sites, sum_dp, n_miss, n_het, sum_frac, sum_ref, sum_alt), sum), .groups = "drop") %>%
  mutate(mean_dp = sum_dp / n_sites, miss = n_miss / n_sites,
         ab = ifelse(n_het >= 20, sum_frac / n_het, NA), ab_reads = sum_ref / (sum_ref + sum_alt))

ph <- read.delim(file.path(root, "input_data/06_Association/PKINGS_ALL_WOB_EXCLUDED.pheno.txt")) %>%
  filter(pheno %in% c(1, 2)) %>% transmute(sample = IID, morph = ifelse(pheno == 2, "BP", "TP"))
md <- read.delim(file.path(root, "input_data/01_Terra/data_model/fish_data_2023.txt"))
md <- md[, c(1, match("POP", names(md)))]   # participant ID and population
names(md) <- c("sample", "POP")

g <- d %>% filter(label == "genome") %>% select(sample, g_dp = mean_dp, g_miss = miss, g_ab = ab)
d <- d %>% left_join(g, "sample") %>% left_join(md, "sample") %>% left_join(ph, "sample") %>%
  mutate(rel_dp = mean_dp / g_dp, d_miss = miss - g_miss, d_ab = ab - g_ab)
write.table(d, file.path(dir, "ref_bias_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## GWAS fish (mixed excluded), BP vs TP, per region class
gw <- d %>% filter(!is.na(morph))
wp <- function(x, y) {  # NA when a group has < 3 usable fish (e.g. ab in small peaks)
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  signif(wilcox.test(x, y, exact = FALSE)$p.value, 3)
}
tab <- gw %>% group_by(label) %>% summarise(
  n_BP = sum(morph == "BP"), n_TP = sum(morph == "TP"),
  sites = first(n_sites),
  dp_BP = mean(mean_dp[morph == "BP"]), dp_TP = mean(mean_dp[morph == "TP"]),
  rel_dp_BP = mean(rel_dp[morph == "BP"]), rel_dp_TP = mean(rel_dp[morph == "TP"]), p_rel_dp = wp(rel_dp[morph == "BP"], rel_dp[morph == "TP"]),
  miss_BP = mean(miss[morph == "BP"]), miss_TP = mean(miss[morph == "TP"]), p_miss = wp(miss[morph == "BP"], miss[morph == "TP"]),
  hets_BP = sum(n_het[morph == "BP"]), hets_TP = sum(n_het[morph == "TP"]),
  ab_BP = mean(ab[morph == "BP"], na.rm = TRUE), ab_TP = mean(ab[morph == "TP"], na.rm = TRUE), p_ab = wp(ab[morph == "BP"], ab[morph == "TP"]),
  d_ab_BP = mean(d_ab[morph == "BP"], na.rm = TRUE), d_ab_TP = mean(d_ab[morph == "TP"], na.rm = TRUE), p_d_ab = wp(d_ab[morph == "BP"], d_ab[morph == "TP"]),
  .groups = "drop") %>%
  mutate(across(where(is.double) & !starts_with("p_"), ~ round(.x, 4)))
write.table(tab, file.path(dir, "ref_bias_BP_vs_TP.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## the three geographically matched comparisons
cmp <- list(BP1_vs_TP1 = list(BP = "BAM", TP = c("APA", "BEN")),
            BP2_vs_TP2 = list(BP = c("BIK", "BAVA"), TP = c("DOG", "DOV")),
            BP3_vs_TP2 = list(BP = "BIR", TP = c("DOG", "DOV")))
by_cmp <- bind_rows(lapply(names(cmp), function(nm) {
  gw %>% filter((morph == "BP" & POP %in% cmp[[nm]]$BP) | (morph == "TP" & POP %in% cmp[[nm]]$TP)) %>%
    group_by(label) %>% summarise(
      comparison = nm, n_BP = sum(morph == "BP"), n_TP = sum(morph == "TP"),
      rel_dp_BP = mean(rel_dp[morph == "BP"]), rel_dp_TP = mean(rel_dp[morph == "TP"]), p_rel_dp = wp(rel_dp[morph == "BP"], rel_dp[morph == "TP"]),
      miss_BP = mean(miss[morph == "BP"]), miss_TP = mean(miss[morph == "TP"]), p_miss = wp(miss[morph == "BP"], miss[morph == "TP"]),
      ab_BP = mean(ab[morph == "BP"], na.rm = TRUE), ab_TP = mean(ab[morph == "TP"], na.rm = TRUE), p_ab = wp(ab[morph == "BP"], ab[morph == "TP"]),
      .groups = "drop")
})) %>% mutate(across(where(is.double) & !starts_with("p_"), ~ round(.x, 4)))
write.table(by_cmp, file.path(dir, "ref_bias_by_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## Per peak, which morph carries the non-reference haplotype? The reference (a Bambomo BP fish) carries the BP
## allele at the Bambomo-origin peaks and the TP allele at the others (orientation key from
## code/14_Wobble_Analysis). Reference bias would penalise the non-reference morph: lower relative depth, more
## missing calls, ref-skewed hets. "mixed" = mean ALT dosage differs by < 1 between BP and TP reference fish.
ori <- read.delim(file.path(root, "output_data/14_Wobble_Analysis/PKINGS_WOB_MULTILOCUS/peak_orientation_key.txt")) %>%
  transmute(label = paste0("peak", peak),
            ref_carries = ifelse(abs(mean_alt_bp_ref - mean_alt_tp_ref) < 1, "mixed", ifelse(bp_allele == "REF", "BP", "TP")))
by_ref <- by_cmp %>% inner_join(ori, "label") %>%
  mutate(nonref = ifelse(ref_carries == "BP", "TP", "BP"),
         sgn = ifelse(nonref == "TP", 1, -1),  # non-reference minus reference-matching morph
         d_rel_dp = round(sgn * (rel_dp_TP - rel_dp_BP), 4), d_miss = round(sgn * (miss_TP - miss_BP), 4),
         d_ab = round(sgn * (ab_TP - ab_BP), 4)) %>%
  select(label, comparison, ref_carries, nonref, d_rel_dp, p_rel_dp, d_miss, p_miss, d_ab, p_ab) %>%
  arrange(label, comparison)
write.table(by_ref, file.path(dir, "ref_bias_peaks_by_ref_haplotype.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## all fish by population (context: other species / distant populations)
by_pop <- d %>% group_by(label, POP) %>%
  summarise(n = n(), dp = mean(mean_dp), rel_dp = mean(rel_dp), miss = mean(miss), ab = mean(ab, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(where(is.double), ~ round(.x, 4)))
write.table(by_pop, file.path(dir, "ref_bias_by_pop.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

options(width = 250)
print(as.data.frame(tab))
print(as.data.frame(by_cmp %>% filter(label %in% c("genome", "hidiv", "peaks_all"))))
print(as.data.frame(by_ref))
