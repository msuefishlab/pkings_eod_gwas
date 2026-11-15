#!/usr/bin/env Rscript
# ============================================================================
# 17_peak_hap_match_audit.R
#
# How BP-like are the "BP-class" haplotypes found OUTSIDE a peak's home
# populations (Bambomo for Peaks 1/4/5; Bikagala/Bavavela/Biroundou for 2/3/6)?
#
# 01_call_peak_haplotypes.R labels each phased haplotype by the nearer of two
# consensus haplotypes (BP vs TP), with no margin, so a haplotype matching the
# BP consensus at 23 of 43 diagnostic SNPs is "BP" just like one matching 43 of
# 43. This script grades those calls by the fraction of diagnostic SNPs that
# carry the BP allele:
#   full           all diagnostic SNPs BP
#   near-complete  >= 90%
#   partial        < 90% (but still nearer the BP consensus)
# The calls themselves are unchanged; the script stops if it does not reproduce
# peak_haplotype_states.csv or the match fractions of per_haplotype.tsv.
#
# Inputs (PKINGS_WOB_MULTILOCUS/, all from 01 and 16):
#   PKINGS_WOB_peak_snps.gt_phased.txt       phased GT, 214 focal fish
#   peak_haplotype_allele_key.csv            BP/TP consensus per peak
#   peak_haplotype_states.csv                per-fish calls (check)
#   nonfocal_haplotypes/per_haplotype.tsv    per-haplotype match, raw-call QC
#   input_data/01_Terra/data_model/fish_data_2024a.txt
#   output_data/08_Peak_Analysis/PKINGS_ALL_WOB_EXCLUDED/..._SNPS_IN_PEAKS_TOP250.vcf.gz
#     (the Fig. S4 SNP set; see README) and the phased VCF in output_data/07_Phasing/;
#     needs bcftools on PATH
#
# Outputs (PKINGS_WOB_MULTILOCUS/peak_hap_zoom/):
#   outside_home_bp_class.tsv   every BP-class haplotype outside its home
#                               populations, with tier, n/N, longest BP run
#   tier_counts.tsv             fish counts by peak x population x tier
#   outside_home_figS4snps.tsv  the same haplotypes scored on the wider Fig. S4
#                               SNP set, with mismatches to the nearest Bambomo
#                               (or BIK/BAVA/BIR) BP-class haplotype
#   peak_hap_zoom_figS4snps.png/.pdf
#                               phased haplotypes at the Fig. S4 SNPs, Peaks 1,
#                               4, 5, coloured BP vs TP allele (not REF/ALT)
#
# Run:  Rscript code/14_Wobble_Analysis/17_peak_hap_match_audit.R
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(patchwork)
})

root <- rprojroot::find_root(rprojroot::has_file(".git") | rprojroot::is_git_root)
setwd(root)

ml_dir <- file.path(root, "output_data", "14_Wobble_Analysis", "PKINGS_WOB_MULTILOCUS")
outdir <- file.path(ml_dir, "peak_hap_zoom")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

NEAR_COMPLETE <- 0.90   # BP-match fraction for "near-complete"
ZOOM_PEAKS    <- c(1, 4, 5)

home_pops <- list(BAM = "BAM", other = c("BIK", "BAVA", "BIR"))
# Peak 6 is a BP2-only locus (Biroundou carries no BP-class copy): its home
# populations, and its BP reference, are Bikagala + Bavavela alone.
in_home_pop <- function(peak, origin, POP)
  (origin == "BAM" & POP %in% home_pops$BAM) |
  (origin == "other" & peak != 6 & POP %in% home_pops$other) |
  (origin == "other" & peak == 6 & POP %in% c("BIK", "BAVA"))

# ---------------------------------------------------------------------------
# 1. Inputs
# ---------------------------------------------------------------------------
fish <- read.table(file.path(root, "input_data", "01_Terra", "data_model", "fish_data_2024a.txt"),
                   header = TRUE, sep = "\t", quote = "") %>%
  transmute(PN = trimws(PN), POP = trimws(POP), Phenotype = trimws(Phenotype), Species = trimws(Species))

key <- read.csv(file.path(ml_dir, "peak_haplotype_allele_key.csv"),
                colClasses = c(bp_consensus = "character", tp_consensus = "character"))
states <- read.csv(file.path(ml_dir, "peak_haplotype_states.csv"))
per_hap <- fread(file.path(ml_dir, "nonfocal_haplotypes", "per_haplotype.tsv"))

gt_path <- file.path(ml_dir, "PKINGS_WOB_peak_snps.gt_phased.txt")
hdr <- strsplit(readLines(gt_path, n = 1), "\t")[[1]]
hdr <- sub(":GT$", "", sub("^#?\\[[0-9]+\\]", "", hdr))
gt <- fread(gt_path, header = FALSE, skip = 1, sep = "\t", colClasses = "character")
setnames(gt, hdr)
gt[, POS := as.integer(trimws(POS))]
gt[, CHROM := trimws(CHROM)]
sample_cols <- setdiff(names(gt), c("CHROM", "POS", "REF", "ALT"))

snps <- read.table(file.path(root, "output_data", "06_Association", "PKINGS_ALL_WOB_EXCLUDED",
                             "PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS.txt"), header = TRUE, sep = "\t")
row_peak <- snps$peak[match(paste(gt$CHROM, gt$POS), paste(snps$chr, snps$ps))]
stopifnot(!any(is.na(row_peak)))

parse_allele <- function(v, which_hap) {
  a <- tstrsplit(v, "|", fixed = TRUE)[[which_hap]]
  a[a == "."] <- NA
  as.integer(a)
}
hapA <- as.matrix(gt[, lapply(.SD, parse_allele, which_hap = 1), .SDcols = sample_cols])
hapB <- as.matrix(gt[, lapply(.SD, parse_allele, which_hap = 2), .SDcols = sample_cols])

# ---------------------------------------------------------------------------
# 2. Per-SNP BP-allele matrix for each peak (diagnostic SNPs only)
# ---------------------------------------------------------------------------
peak_long <- function(pk) {
  k <- key[key$peak == pk, ]
  rows <- which(row_peak == pk)
  bp_cons <- as.integer(strsplit(k$bp_consensus, "")[[1]])
  tp_cons <- as.integer(strsplit(k$tp_consensus, "")[[1]])
  stopifnot(length(bp_cons) == length(rows))
  diag <- which(bp_cons != tp_cons)
  stopifnot(length(diag) == k$n_informative)
  r <- rows[diag]
  bind_rows(
    lapply(c(A = "A", B = "B"), function(h) {
      M <- if (h == "A") hapA else hapB
      is_bp <- M[r, , drop = FALSE] == bp_cons[diag]
      as_tibble(is_bp) %>%
        mutate(snp = seq_along(r), pos = gt$POS[r]) %>%
        pivot_longer(all_of(sample_cols), names_to = "PN", values_to = "bp") %>%
        mutate(hap = h)
    })
  ) %>% mutate(peak = pk, origin = k$origin, n_diag = length(diag))
}
long <- bind_rows(lapply(sort(unique(row_peak)), peak_long)) %>%
  left_join(fish, by = "PN")

hap_sum <- long %>%
  group_by(peak, origin, PN, POP, Phenotype, hap, n_diag) %>%
  summarise(n_bp = sum(bp, na.rm = TRUE), n_tp = sum(!bp, na.rm = TRUE), .groups = "drop") %>%
  mutate(bp_match_frac = n_bp / n_diag,
         hap_class = case_when(n_bp > n_tp ~ "BP", n_bp < n_tp ~ "TP", TRUE ~ "ambig"))

# Check 1: reproduce 01's per-haplotype classes.
chk01 <- states %>%
  select(PN, peak, hapA_class, hapB_class) %>%
  pivot_longer(c(hapA_class, hapB_class), names_to = "hap", values_to = "class01") %>%
  mutate(hap = sub("hap(.)_class", "\\1", hap)) %>%
  inner_join(hap_sum, by = c("PN", "peak", "hap"))
stopifnot(nrow(chk01) == nrow(hap_sum), all(chk01$class01 == chk01$hap_class))

# Check 2: reproduce 16's match fractions for the focal fish.
chk16 <- per_hap %>%
  select(PN, peak, hap, bp_match_frac16 = bp_match_frac) %>%
  inner_join(hap_sum, by = c("PN", "peak", "hap"))
stopifnot(nrow(chk16) == nrow(hap_sum),
          all(abs(chk16$bp_match_frac16 - chk16$bp_match_frac) < 0.001))
cat(sprintf("Reproduced %d haplotype calls (01) and match fractions (16).\n", nrow(hap_sum)))

# ---------------------------------------------------------------------------
# 3. Grade BP-class haplotypes outside their home populations
# ---------------------------------------------------------------------------
tier_of <- function(f) case_when(f == 1 ~ "full", f >= NEAR_COMPLETE ~ "near-complete", TRUE ~ "partial")

outside <- per_hap %>%
  filter(focal) %>%
  filter(!in_home_pop(peak, origin, POP),
         hap_class == "BP") %>%
  mutate(tier = tier_of(bp_match_frac)) %>%
  select(peak, chr, origin, PN, POP, Phenotype, hap, n_diag, n_bp_match, bp_match_frac, tier,
         longest_bp_run_snps, longest_bp_run_bp, raw_called_frac, raw_phased_concord, median_dp) %>%
  arrange(peak, POP, desc(bp_match_frac), PN)
write_tsv(outside, file.path(outdir, "outside_home_bp_class.tsv"))

# One row per fish: its best BP-class copy decides the tier.
tier_rank <- c("full", "near-complete", "partial")
n_fish <- states %>% count(peak, POP, Phenotype, name = "n_fish")
tiers <- outside %>%
  group_by(peak, POP, Phenotype, PN) %>%
  summarise(tier = tier_rank[min(match(tier, tier_rank))], .groups = "drop") %>%
  count(peak, POP, Phenotype, tier) %>%
  pivot_wider(names_from = tier, values_from = n, values_fill = 0) %>%
  left_join(n_fish, by = c("peak", "POP", "Phenotype")) %>%
  select(peak, POP, Phenotype, n_fish, any_of(tier_rank))
write_tsv(tiers, file.path(outdir, "tier_counts.tsv"))

cat("\n=== BP-class haplotypes outside home populations, fish counts by best copy ===\n")
print(as.data.frame(tiers))

# ---------------------------------------------------------------------------
# 4. Fig. S4 SNP set: phased haplotypes at the SNPs Fig. S4 plots
#    (08_Peak_Analysis TOP250 VCF: the extended -log10 p > 3 peak SNPs, wider
#    than the -log10 p > 5 peak SNPs used above).
#    Same reference pools and consensus rule as 01, but diagnostic SNPs are
#    re-derived on this set. Each haplotype outside Bambomo is also compared
#    with every Bambomo BP-class haplotype directly (nearest-neighbour
#    mismatches over all window SNPs), so a carrier of a minority Bambomo
#    subtype is not penalised for differing from the Bambomo majority.
# ---------------------------------------------------------------------------
figs4_vcf <- file.path(root, "output_data", "08_Peak_Analysis", "PKINGS_ALL_WOB_EXCLUDED",
                       "PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS_TOP250.vcf.gz")
phased_vcf <- file.path(root, "output_data", "07_Phasing",
                        "MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_SNPS_ONLY.phased.renamed.vcf.gz")
sites_f <- tempfile(fileext = ".tsv"); s4_gt_f <- tempfile(fileext = ".tsv")
system2("bcftools", c("query", "-f", "'%CHROM\\t%POS\\n'", figs4_vcf), stdout = sites_f)
system(sprintf("bcftools view -T %s -R %s %s 2>/dev/null | bcftools query -H -f '%%CHROM\\t%%POS[\\t%%GT]\\n' > %s",
               sites_f, sites_f, phased_vcf, s4_gt_f))
hdr4 <- sub(":GT$", "", sub("^#?\\s*\\[[0-9]+\\]", "", strsplit(readLines(s4_gt_f, n = 1), "\t")[[1]]))
g4 <- fread(s4_gt_f, header = FALSE, skip = 1, sep = "\t", colClasses = "character")
setnames(g4, hdr4)
g4[, POS := as.integer(POS)]
n_sites <- nrow(fread(sites_f, header = FALSE))
stopifnot(nrow(g4) == n_sites)   # every Fig. S4 SNP is in the phased VCF
s4_samples <- intersect(setdiff(names(g4), c("CHROM", "POS")), fish$PN)
A4 <- as.matrix(g4[, lapply(.SD, parse_allele, which_hap = 1), .SDcols = s4_samples])
B4 <- as.matrix(g4[, lapply(.SD, parse_allele, which_hap = 2), .SDcols = s4_samples])
bed <- read.table(file.path(root, "output_data", "06_Association", "PKINGS_ALL_WOB_EXCLUDED",
                            "PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed"),
                  col.names = c("chr_num", "start", "end", "peak", "maxp", "maxp_pos"))
g4_peak <- bed$peak[match(g4$CHROM, paste0("chr", bed$chr_num))]
stopifnot(!any(is.na(g4_peak)))

bam_bp_ids   <- intersect(fish$PN[fish$POP == "BAM" & fish$Phenotype == "BP"], s4_samples)
other_bp_ids <- intersect(fish$PN[fish$POP %in% home_pops$other & fish$Phenotype == "BP"], s4_samples)
tp_ref_ids   <- intersect(fish$PN[fish$POP %in% c("APA", "BEN", "DOV", "DOG") & fish$Phenotype == "TP"], s4_samples)

# A non-home haplotype "is a Bambomo BP haplotype" when it lies within
# NN_FRAC of the window's SNPs (at least 1) of some Bambomo BP-class haplotype.
NN_FRAC <- 0.02

src_ids <- list(BAM = bam_bp_ids,
                BIKBAVA = intersect(fish$PN[fish$POP %in% c("BIK", "BAVA") & fish$Phenotype == "BP"], s4_samples),
                BIR = intersect(fish$PN[fish$POP == "BIR" & fish$Phenotype == "BP"], s4_samples))
bp_ref_for <- function(pk, origin) {
  if (origin == "BAM") bam_bp_ids else if (pk == 6) src_ids$BIKBAVA else other_bp_ids
}

window_hap <- function(pk) {
  r <- which(g4_peak == pk)
  H <- cbind(A4[r, , drop = FALSE], B4[r, , drop = FALSE])
  colnames(H) <- c(paste0(s4_samples, "|A"), paste0(s4_samples, "|B"))
  H
}
ref_cols <- function(ids) c(paste0(ids, "|A"), paste0(ids, "|B"))

figs4_peak <- function(pk) {
  r <- which(g4_peak == pk)
  origin <- key$origin[key$peak == pk]
  bp_ref <- bp_ref_for(pk, origin)
  H <- window_hap(pk)
  bp_freq <- rowMeans(H[, ref_cols(bp_ref)], na.rm = TRUE)
  tp_freq <- rowMeans(H[, ref_cols(tp_ref_ids)], na.rm = TRUE)
  bp_cons <- as.integer(bp_freq >= 0.5)
  tp_cons <- as.integer(tp_freq >= 0.5)
  diag <- which(bp_cons != tp_cons)
  is_bp <- H[diag, , drop = FALSE] == bp_cons[diag]
  n_bp <- colSums(is_bp, na.rm = TRUE); n_tp <- colSums(!is_bp, na.rm = TRUE)
  # BP allele at the index SNP (01 orientation rule: more common in the BP pool).
  idx <- which(g4$POS[r] == bed$maxp_pos[bed$peak == pk])
  stopifnot(length(idx) == 1)
  idx_bp_allele <- as.integer(bp_freq[idx] > tp_freq[idx])
  hs <- tibble(id = colnames(H), PN = sub("\\|.$", "", id), hap = sub("^.*\\|", "", id),
               peak = pk, origin = origin, n_snps = length(r), n_diag = length(diag),
               n_bp = n_bp, bp_match_frac = n_bp / length(diag),
               hap_class = case_when(n_bp > n_tp ~ "BP", n_bp < n_tp ~ "TP", TRUE ~ "ambig"),
               index_bp = as.integer(H[idx, ] == idx_bp_allele)) %>%
    left_join(fish, by = "PN")
  k <- max(1, round(NN_FRAC * length(r)))
  # Nearest BP-class haplotype of a reference set, over all window SNPs (self excluded).
  nn_to <- function(ref_ids) {
    ref <- intersect(ref_cols(ref_ids), hs$id[hs$hap_class == "BP"])
    t(vapply(hs$id, function(h) {
      rr <- setdiff(ref, h)
      if (length(rr) == 0) return(c(NA_real_, NA_real_, 0))
      d <- colSums(H[, rr, drop = FALSE] != H[, h], na.rm = TRUE)
      c(min(d), match(names(which.min(d)), colnames(H)), sum(d <= k))
    }, numeric(3)))
  }
  home_ids <- if (origin == "BAM") src_ids$BAM else if (pk == 6) src_ids$BIKBAVA else c(src_ids$BIKBAVA, src_ids$BIR)
  nn <- nn_to(home_ids)
  out <- hs %>% mutate(nn_home_dist = nn[, 1], nn_home_hap = colnames(H)[nn[, 2]],
                       n_home_within_k = nn[, 3], k = k,
                       home_type = hap_class == "BP" & !is.na(nn_home_dist) & nn_home_dist <= k)
  for (src in names(src_ids)) out[[paste0("nn_", src)]] <- nn_to(src_ids[[src]])[, 1]
  out
}
figs4 <- bind_rows(lapply(sort(unique(g4_peak)), figs4_peak))

# Per-SNP BP-allele calls (diagnostic SNPs, origin-specific consensus) for the figures.
figs4_long <- bind_rows(lapply(sort(unique(g4_peak)), function(pk) {
  r <- which(g4_peak == pk)
  H <- window_hap(pk)
  bp_ref <- bp_ref_for(pk, key$origin[key$peak == pk])
  cons <- function(ids) as.integer(rowMeans(H[, ref_cols(ids)], na.rm = TRUE) >= 0.5)
  bp_cons <- cons(bp_ref); diag <- which(bp_cons != cons(tp_ref_ids))
  as_tibble(H[diag, , drop = FALSE] == bp_cons[diag]) %>%
    mutate(snp = seq_along(diag), pos = g4$POS[r][diag], peak = pk) %>%
    pivot_longer(-c(snp, pos, peak), names_to = "id", values_to = "bp")
}))

is_home <- function(peak, origin, POP, Phenotype) Phenotype == "BP" & in_home_pop(peak, origin, POP)

outside4 <- figs4 %>%
  filter(!in_home_pop(peak, origin, POP), hap_class == "BP") %>%
  select(peak, PN, POP, Species, Phenotype, hap, n_snps, n_diag, n_bp, bp_match_frac, index_bp,
         nn_home_dist, nn_home_hap, n_home_within_k, k, home_type, nn_BAM, nn_BIKBAVA, nn_BIR) %>%
  arrange(peak, POP, nn_home_dist, desc(bp_match_frac))
write_tsv(outside4, file.path(outdir, "outside_home_figS4snps.tsv"))

# Population summary: every population and species, every peak.
sharing <- figs4 %>%
  group_by(peak, origin, POP, Species, Phenotype) %>%
  summarise(n_fish = n_distinct(PN), n_snps = first(n_snps), k = first(k),
            index_bp_haps = sum(index_bp, na.rm = TRUE), n_haps = n(),
            bp_class_haps = sum(hap_class == "BP"),
            home_type_haps = sum(home_type), home_type_fish = n_distinct(PN[home_type]),
            min_nn_home = suppressWarnings(min(nn_home_dist[hap_class == "BP"])),
            bir_within_k_of_bikbava = sum(hap_class == "BP" & nn_BIKBAVA <= k, na.rm = TRUE),
            bikbava_within_k_of_bir = sum(hap_class == "BP" & nn_BIR <= k, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(home = is_home(peak, origin, POP, Phenotype), home_pop = in_home_pop(peak, origin, POP),
         min_nn_home = ifelse(is.infinite(min_nn_home), NA, min_nn_home)) %>%
  arrange(peak, desc(home), POP, Phenotype)
write_tsv(sharing, file.path(outdir, "sharing_by_population_figS4snps.tsv"))

cat("\n=== Fig. S4 SNP set: home-type haplotypes outside the home populations ===\n")
print(sharing %>% filter(!home_pop, bp_class_haps > 0 | index_bp_haps > 0) %>%
        select(peak, POP, Phenotype, n_fish, k, index_bp_haps, bp_class_haps, home_type_haps,
               home_type_fish, min_nn_home) %>% as.data.frame())
cat("\n=== Calibration: home BP-class haplotypes, mismatches to the nearest other home BP-class haplotype ===\n")
print(figs4 %>% filter(is_home(peak, origin, POP, Phenotype), hap_class == "BP") %>%
        group_by(peak, n_snps, k) %>%
        summarise(n_haps = n(), median = median(nn_home_dist), q90 = quantile(nn_home_dist, 0.9),
                  max = max(nn_home_dist), .groups = "drop") %>% as.data.frame())
cat("\n=== BP2 vs BP3 at the other-origin peaks (BP-class haps within k of the other group) ===\n")
print(sharing %>% filter(origin == "other", home) %>%
        select(peak, POP, n_fish, k, bp_class_haps, bir_within_k_of_bikbava, bikbava_within_k_of_bir) %>%
        as.data.frame())
print(figs4 %>% filter(origin == "other", POP == "BIR", Phenotype == "BP", hap_class == "BP") %>%
        group_by(peak) %>%
        summarise(bir_bp_haps = n(), nn_BIKBAVA_min = min(nn_BIKBAVA), nn_BIKBAVA_median = median(nn_BIKBAVA),
                  nn_BIKBAVA_max = max(nn_BIKBAVA), .groups = "drop") %>% as.data.frame())

# ---------------------------------------------------------------------------
# 5. Zoom figures on the Fig. S4 SNP set: Peaks 1/4/5 and Peaks 2/3/6
# ---------------------------------------------------------------------------
grp_of <- function(POP, Ph) case_when(
  POP == "BAM" & Ph == "BP"                ~ "Bambomo BP",
  POP %in% c("BIK", "BAVA") & Ph == "BP"   ~ "Bikagala + Bavavela BP",
  POP == "BIR" & Ph == "BP"                ~ "Biroundou BP",
  POP == "APA" & Ph == "TP"                ~ "Apassa TP",
  POP == "APA" & Ph == "Wobble"            ~ "Apassa mixed",
  POP == "BEN" & Ph == "TP"                ~ "Bengue TP",
  POP %in% c("DOV", "DOG")                 ~ "Douvalou + Douengi TP",
  POP %in% c("MOV", "IVI")                 ~ "Mouvanga + Ivindo TP",
  POP == "COB"                             ~ "Cobia BP + TP",
  POP %in% c("PCUR", "PNTO", "PSN2")       ~ "P. curvifrons, ntotom, SN2 (BP)",
  POP %in% c("PLON", "PPAR", "PGAB", "PHOP") ~ "P. longicaudatus, parsegnis,\ngabonensis, hopkinsi (BP)",
  POP == "PMAG"                            ~ "P. magnostipes",
  POP %in% c("PSZA", "PVAD")               ~ "P. szaboi, vadamans (TP)")
grp_rest <- c("Apassa TP", "Apassa mixed", "Bengue TP", "Douvalou + Douengi TP", "Mouvanga + Ivindo TP",
              "Cobia BP + TP", "P. curvifrons, ntotom, SN2 (BP)",
              "P. longicaudatus, parsegnis,\ngabonensis, hopkinsi (BP)", "P. magnostipes",
              "P. szaboi, vadamans (TP)")
fig_sets <- list(
  BAM = list(peaks = c(1, 4, 5), home_label = "Bambomo",
             levels = c("Bambomo BP", "Biroundou BP", "Bikagala + Bavavela BP", grp_rest),
             file = "peak_hap_zoom_figS4snps",
             title = "Phased haplotypes at the Fig. S4 SNPs, Bambomo-origin peaks (1, 4, 5)"),
  other = list(peaks = c(2, 3, 6), home_label = "BIK/BAVA/BIR",
               levels = c("Bikagala + Bavavela BP", "Biroundou BP", "Bambomo BP", grp_rest),
               file = "peak_hap_zoom_figS4snps_peaks236",
               title = "Phased haplotypes at the Fig. S4 SNPs, Bikagala/Bavavela/Biroundou-origin peaks (2, 3, 6)"))
MIN_ROWS <- 26   # small groups are padded with blank rows so their notes fit

col_bp <- "#1F4E9C"   # BP allele
col_tp <- "#E3E3E3"   # TP allele

name_or_count <- function(PN, POP, dist) {
  if (n_distinct(PN) <= 4) return(paste(sprintf("%s (%d)", PN, dist), collapse = "\n  "))
  by_pop <- tibble(PN, POP) %>% distinct() %>% count(POP) %>% mutate(t = paste(POP, n))
  sprintf("%d haps in %d fish\n  (%s)", length(PN), n_distinct(PN), paste(by_pop$t, collapse = ", "))
}

# Per-group note. Non-home groups: home-type haplotypes (named, with mismatches
# to the nearest home BP-class haplotype) and other BP-like haplotypes
# (counted). Biroundou at Peaks 2/3/6: its BP-class haplotypes within k of a
# Bikagala/Bavavela one.
group_note <- function(hs, home_label) {
  kk <- first(hs$k); snps_lab <- first(hs$n_snps)
  x <- hs %>% filter(hap_class == "BP", !home_pop)
  home_t <- x %>% filter(home_type) %>% group_by(group) %>%
    summarise(line = sprintf("%s BP type (<=%d/%d):\n  %s", home_label, kk, snps_lab,
                             name_or_count(PN, POP, nn_home_dist)), ord = 1, .groups = "drop")
  oth <- x %>% filter(!home_type) %>% group_by(group) %>%
    summarise(line = sprintf("BP-like, not %s:\n  %d hap%s in %d fish (%s mismatches)", home_label,
                             n(), if (n() == 1) "" else "s", n_distinct(PN),
                             if (min(nn_home_dist) == max(nn_home_dist)) min(nn_home_dist)
                             else paste0(min(nn_home_dist), "-", max(nn_home_dist))),
              ord = 2, .groups = "drop")
  bir <- hs %>% filter(home, POP == "BIR", hap_class == "BP") %>% group_by(group) %>%
    summarise(line = sprintf("BP-class haps within %d/%d of a\n  Bikagala/Bavavela BP hap: %d of %d",
                             kk, snps_lab, sum(nn_BIKBAVA <= kk, na.rm = TRUE), n()), ord = 1, .groups = "drop")
  shown <- levels(droplevels(hs$group))
  none <- tibble(group = setdiff(shown[!shown %in% hs$group[hs$home]], as.character(x$group))) %>%
    mutate(group = factor(group, levels = levels(hs$group)), line = "no BP-like haplotype", ord = 3)
  bind_rows(home_t, oth, bir, none) %>% arrange(group, ord) %>% group_by(group) %>%
    summarise(note = paste(line, collapse = "\n"), .groups = "drop")
}

zoom_panel <- function(pk, fs, show_strips) {
  hs <- figs4 %>% filter(peak == pk) %>%
    mutate(group = factor(grp_of(POP, Phenotype), levels = fs$levels),
           home = is_home(peak, origin, POP, Phenotype), home_pop = in_home_pop(peak, origin, POP)) %>%
    filter(!is.na(group)) %>%
    arrange(group, if_else(home, 0, nn_home_dist), desc(bp_match_frac), PN, hap)
  # Pad small groups with blank rows (bottom of the group).
  pad <- hs %>% count(group) %>% filter(n < MIN_ROWS) %>%
    rowwise() %>% reframe(group = group, id = sprintf("pad|%s|%02d", group, seq_len(MIN_ROWS - n)))
  ord <- bind_rows(hs %>% select(group, id), pad) %>% arrange(group)
  ord <- ord %>% mutate(row = factor(id, levels = rev(id)))
  hs <- hs %>% left_join(ord %>% select(id, row), by = "id")
  d <- figs4_long %>% filter(peak == pk) %>% inner_join(hs %>% select(id, row, group), by = "id")
  d_pad <- tidyr::crossing(ord %>% filter(startsWith(id, "pad|")), snp = range(d$snp)) %>%
    mutate(bp = NA)
  notes <- group_note(hs, fs$home_label) %>%
    left_join(ord %>% group_by(group) %>% summarise(row = first(row), .groups = "drop"), by = "group") %>%
    mutate(x = max(d$snp) * 1.035 + 1.5)
  ticks <- hs %>% filter(home_type, !home_pop) %>%
    mutate(x = max(d$snp) + 0.5 + max(d$snp) * 0.0125)
  span_kb <- (max(d$pos) - min(d$pos)) / 1e3
  cur <- bed[bed$peak == pk, ]
  cur_snps <- range(d$snp[d$pos >= cur$start & d$pos <= cur$end])
  ggplot(d, aes(snp, row, fill = bp)) +
    geom_raster() +
    geom_blank(data = d_pad) +
    annotate("segment", x = cur_snps[1] - 0.5, xend = cur_snps[2] + 0.5, y = Inf, yend = Inf,
             colour = "#D55E00", linewidth = 1.6) +
    geom_tile(data = ticks, aes(x = x, y = row), inherit.aes = FALSE,
              width = max(d$snp) * 0.02, height = 1, fill = "black") +
    geom_text(data = notes, aes(x = x, y = row, label = note), inherit.aes = FALSE,
              hjust = 0, vjust = 1, size = 3.3, lineheight = 0.9, colour = "black") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = c(`TRUE` = col_bp, `FALSE` = col_tp),
                      labels = c(`TRUE` = "BP allele", `FALSE` = "TP allele"),
                      name = NULL, na.value = "white", na.translate = FALSE) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = sprintf("Peak %d  (%s; %d/%d SNPs diagnostic, %.0f kb)",
                         pk, key$chr[key$peak == pk], first(hs$n_diag), first(hs$n_snps), span_kb),
         x = "Diagnostic SNPs, in genome order", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(),
          panel.spacing.y = unit(3, "pt"),
          strip.text.y.left = if (show_strips) element_text(angle = 0, hjust = 1, size = 10.5) else element_blank(),
          strip.placement = "outside",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 11.5, face = "bold"),
          plot.margin = margin(4, 215, 4, 4),
          legend.position = "bottom")
}

for (fs in fig_sets) {
  ref <- fs$levels[1]
  p <- wrap_plots(lapply(seq_along(fs$peaks), function(i) zoom_panel(fs$peaks[i], fs, i == 1)), nrow = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = fs$title,
      subtitle = paste0(
        "BP/TP allele = ", fs$home_label, " BP vs Apassa/Bengue/Douvalou/Douengi TP consensus. ",
        "Orange bar: SNPs inside the current peak interval. kb = span of the diagnostic SNPs.\n",
        "Rows sorted by mismatches to the nearest ", fs$home_label, " BP-class haplotype over all window SNPs. ",
        "Black tick, ", fs$home_label, " BP type: within 2% of the window's SNPs of one (<= k/N in the note);\n",
        "(n) = mismatches to it. White = padding for small groups. Bambomo mixed/TP and Bikagala/Bavavela mixed are not shown."),
      theme = theme(plot.title = element_text(size = 14, face = "bold"),
                    plot.subtitle = element_text(size = 10.5),
                    legend.position = "bottom"))
  ggsave(file.path(outdir, paste0(fs$file, ".png")), p, width = 20, height = 17, dpi = 300, bg = "white")
  ggsave(file.path(outdir, paste0(fs$file, ".pdf")), p, width = 20, height = 17, bg = "white")
}
cat(sprintf("\nWrote to %s:\n  %s\n", outdir,
            paste(c("outside_home_bp_class.tsv", "tier_counts.tsv", "outside_home_figS4snps.tsv",
                    "sharing_by_population_figS4snps.tsv",
                    paste0(vapply(fig_sets, `[[`, "", "file"), ".png")), collapse = "\n  ")))

# Per-haplotype scores and per-SNP calls for 18_figS4_replacement.R.
saveRDS(list(figs4 = figs4, figs4_long = figs4_long, bed = bed, key = key),
        file.path(outdir, "figs4_haplotypes.rds"))
