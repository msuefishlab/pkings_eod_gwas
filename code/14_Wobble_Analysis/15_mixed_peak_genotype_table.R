#!/usr/bin/env Rscript
# ============================================================================
# 15_mixed_peak_genotype_table.R
#
# Numbers for the Supplementary Note on "mixed" (wobble) individuals: their genotypes at
# the six GWAS peaks, heterozygosity vs BP fish, ADMIXTURE ancestry, and the
# fish with histology-confirmed electrocyte anatomy.
#
# Inputs (relative to repo root):
#   output_data/14_Wobble_Analysis/PKINGS_WOB_MULTILOCUS/peak_haplotype_states.csv
#     (phased per-peak states from 01_call_peak_haplotypes.R)
#   output_data/13_Structure/kings_only/*.kings_only.{8,10}.Q + .fam (ADMIXTURE)
#   input_data/01_Terra/data_model/fish_data_2024a.txt
#
# Outputs (PKINGS_WOB_MULTILOCUS/):
#   mixed_peak_state_counts.tsv  population x phenotype x peak: TP_hom/het/BP_hom (Table S2)
#   mixed_het_tests.tsv          BAM mixed vs BAM BP heterozygote fraction at peaks 4 and 5
#   mixed_admixture.tsv          per-fish ancestry outside the population's own cluster, K = 8 and 10
#   anatomy_fish_states.tsv      peak states of the 7 histology-confirmed fish
#
# Run:  Rscript code/14_Wobble_Analysis/15_mixed_peak_genotype_table.R
# ============================================================================

suppressPackageStartupMessages(library(tidyverse))

root <- rprojroot::find_root(rprojroot::is_git_root)
setwd(root)

outdir <- file.path("output_data", "14_Wobble_Analysis", "PKINGS_WOB_MULTILOCUS")
states <- read_csv(file.path(outdir, "peak_haplotype_states.csv"), show_col_types = FALSE)
pops <- c("BAM", "APA", "BEN", "BIK", "BAVA", "BIR", "DOG", "DOV")

# ---- Table S2: peak states by population x phenotype ----------------------
counts <- states %>%
  filter(POP %in% pops) %>%
  count(POP, Phenotype, peak, origin, peak_state) %>%
  pivot_wider(names_from = peak_state, values_from = n, values_fill = 0) %>%
  left_join(states %>% filter(POP %in% pops) %>% distinct(PN, POP, Phenotype) %>% count(POP, Phenotype, name = "n_fish"),
            by = c("POP", "Phenotype")) %>%
  mutate(POP = factor(POP, levels = pops),
         Phenotype = factor(Phenotype, levels = c("BP", "Wobble", "TP"))) %>%
  arrange(POP, Phenotype, peak) %>%
  select(POP, Phenotype, n_fish, peak, origin, any_of(c("TP_hom", "het", "BP_hom", "ambiguous")))
write_tsv(counts, file.path(outdir, "mixed_peak_state_counts.tsv"))

# Full BP genotype at the three Bambomo-origin peaks (BP haplotype present at 1, 4 and 5).
wide <- states %>%
  select(PN, POP, Phenotype, peak, peak_state) %>%
  pivot_wider(names_from = peak, values_from = peak_state, names_prefix = "P")
bam_full <- wide %>%
  filter(POP == "BAM") %>%
  mutate(full_bp = P1 != "TP_hom" & P4 != "TP_hom" & P5 != "TP_hom")
cat("BAM fish carrying a BP haplotype at peaks 1, 4 and 5:\n")
print(count(bam_full, Phenotype, full_bp))
cat("BAM mixed fish without the full BP genotype:\n")
print(filter(bam_full, Phenotype == "Wobble", !full_bp))

# ---- heterozygosity, BAM mixed vs BAM BP ----------------------------------
het_tests <- map_dfr(c(4, 5), function(p) {
  d <- states %>% filter(POP == "BAM", Phenotype %in% c("BP", "Wobble"), peak == p)
  tab <- table(factor(d$Phenotype, levels = c("BP", "Wobble")), d$peak_state == "het")
  tibble(peak = p,
         het_bp = tab["BP", "TRUE"], n_bp = sum(tab["BP", ]),
         het_mixed = tab["Wobble", "TRUE"], n_mixed = sum(tab["Wobble", ]),
         fisher_p = fisher.test(tab)$p.value)
})
write_tsv(het_tests, file.path(outdir, "mixed_het_tests.tsv"))
cat("\nHeterozygote fraction, BAM BP vs BAM mixed:\n")
print(het_tests)

# ---- ADMIXTURE: ancestry from populations of the opposite phenotype -------
# Each cluster is assigned to the population with the highest mean Q among
# non-mixed fish, and takes that population's majority phenotype (BP or TP).
# opposite_ancestry = summed Q over clusters whose phenotype differs from the
# majority phenotype of the fish's own population. A TP x BP F1 or F2 would
# score ~0.5 and a first backcross ~0.25. (1 - own-cluster Q is not used: at
# K = 10 Apassa itself splits between an Apassa and a Bengue-like cluster.)
meta <- read_tsv("input_data/01_Terra/data_model/fish_data_2024a.txt", show_col_types = FALSE,
                 col_types = cols(.default = "c")) %>%
  transmute(PN = trimws(PN), POP = trimws(POP), Phenotype = trimws(Phenotype))
pop_pheno <- meta %>% filter(Phenotype %in% c("BP", "TP")) %>%
  count(POP, Phenotype) %>% group_by(POP) %>% slice_max(n, n = 1, with_ties = FALSE) %>%
  select(POP, pop_pheno = Phenotype)
adm_dir <- file.path("output_data", "13_Structure", "kings_only")
fam <- read.table(list.files(adm_dir, pattern = "kings_only\\.fam$", full.names = TRUE))$V2
admix <- map_dfr(c(8, 10), function(k) {
  q <- read.table(list.files(adm_dir, pattern = sprintf("kings_only\\.%d\\.Q$", k), full.names = TRUE))
  stopifnot(nrow(q) == length(fam))
  long <- as_tibble(q) %>% mutate(PN = fam) %>% left_join(meta, by = "PN") %>%
    left_join(pop_pheno, by = "POP") %>%
    pivot_longer(starts_with("V"), names_to = "cluster", values_to = "Q")
  cl <- long %>% filter(Phenotype != "Wobble") %>%
    group_by(cluster, POP, pop_pheno) %>% summarise(mQ = mean(Q), .groups = "drop") %>%
    group_by(cluster) %>% slice_max(mQ, n = 1, with_ties = FALSE) %>%
    select(cluster, cluster_pop = POP, cluster_pheno = pop_pheno)
  long %>% left_join(cl, by = "cluster") %>%
    group_by(PN, POP, Phenotype, pop_pheno) %>%
    summarise(K = k,
              opposite_ancestry = sum(Q[cluster_pheno != pop_pheno]),
              bam_cluster = sum(Q[cluster_pop == "BAM"]),
              minor = 1 - max(Q), .groups = "drop")
})
write_tsv(admix, file.path(outdir, "mixed_admixture.tsv"))
cat("\nAncestry from opposite-phenotype populations, by population x phenotype:\n")
admix %>% filter(POP %in% c("BAM", "APA", "BIK", "BAVA")) %>%
  group_by(K, POP, Phenotype) %>%
  summarise(n = n(), median_opp = median(opposite_ancestry), max_opp = max(opposite_ancestry),
            median_bam = median(bam_cluster), max_bam = max(bam_cluster), .groups = "drop") %>%
  print(n = Inf)
cat("\nMax opposite-phenotype ancestry over all mixed fish:\n")
admix %>% filter(Phenotype == "Wobble") %>% group_by(K) %>%
  summarise(n = n(), max = max(opposite_ancestry), who = PN[which.max(opposite_ancestry)]) %>% print()
cat("\nApassa mixed vs Apassa TP, Bambomo-cluster ancestry (Wilcoxon):\n")
admix %>% filter(POP == "APA") %>% group_by(K) %>%
  summarise(p = wilcox.test(bam_cluster[Phenotype == "Wobble"], bam_cluster[Phenotype == "TP"], exact = FALSE)$p.value,
            med_mixed = median(bam_cluster[Phenotype == "Wobble"]), med_tp = median(bam_cluster[Phenotype == "TP"])) %>%
  print()
cat("\nCV error by K:\n")
list.files(adm_dir, pattern = "kings_only\\.[0-9]+\\.out$", full.names = TRUE) %>%
  map_chr(~ grep("CV error", readLines(.x), value = TRUE)) %>% sort() %>% cat(sep = "\n")

# ---- histology-confirmed fish ---------------------------------------------
# Anatomy calls as in code/02_Phenotyping/p0_statistical_detection.Rmd (chunk `anatomy`).
anat <- tribble(
  ~PN,        ~anatomy,
  "APA_6737", "TP", "APA_6675", "TP", "APA_6676", "TP", "BAM_6597", "TP",
  "BAM_6494", "Mosaic", "APA_6681", "Mosaic",
  "BAM_6497", "BP"
)
anat_states <- anat %>% left_join(wide, by = "PN")
write_tsv(anat_states, file.path(outdir, "anatomy_fish_states.tsv"))
cat("\nHistology-confirmed fish (EOD label = Phenotype):\n")
print(anat_states)
