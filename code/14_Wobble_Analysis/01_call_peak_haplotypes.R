#!/usr/bin/env Rscript
# ============================================================================
# 01_call_peak_haplotypes.R
#
# Post-GWAS multilocus HAPLOTYPE calling at the six BP-vs-TP GWAS peaks.
#
# For each individual, uses the PHASED genotypes at each peak's SNPs to assign
# each of the two chromosome copies to a BP-like or TP-like haplotype allele,
# then collapses to a per-peak diploid state (BP_hom / het / TP_hom / ambiguous).
# Within-peak phasing collapses many correlated SNPs into one BP-vs-TP haplotype
# call per copy (rather than the dosage at the index SNP alone).
#
# The six peaks sit on six different chromosomes, so there is no cross-peak
# phase. "Recombinant across peaks" therefore means a noncanonical COMBINATION
# of independent per-peak states, which downstream scripts test.
#
# Inputs (all relative to repo root):
#   output_data/14_Wobble_Analysis/PKINGS_WOB_MULTILOCUS/PKINGS_WOB_peak_snps.gt_phased.txt
#   output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed
#   output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS.txt
#   input_data/01_Terra/data_model/fish_data_2024a.txt
#
# Outputs (to PKINGS_WOB_MULTILOCUS/):
#   peak_haplotype_states.csv         long: PN, peak, chr, origin, bp_hap_count, peak_state, ...
#   peak_haplotype_allele_key.csv     per-peak haplotype-allele definition summary
#   peak_haplotype_concordance.csv    agreement vs index-SNP dosage state (robustness)
#
# Run:  Rscript code/14_Wobble_Analysis/01_call_peak_haplotypes.R
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

root <- rprojroot::find_root(rprojroot::has_file(".git") | rprojroot::is_git_root)
setwd(root)

outdir <- file.path(root, "output_data", "14_Wobble_Analysis", "PKINGS_WOB_MULTILOCUS")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- tunables -------------------------------------------------------------
AMBIG_MARGIN <- 0   # a haplotype is 'ambiguous' if |d_bp - d_tp| <= this many
                    # informative SNPs (0 = only exact ties are ambiguous).

# ---------------------------------------------------------------------------
# 1. Peak definitions + peak->origin map
# ---------------------------------------------------------------------------
peaks_path <- file.path(root, "output_data", "06_Association",
                        "PKINGS_ALL_WOB_EXCLUDED",
                        "PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed")
df.bed <- read.table(peaks_path, header = FALSE)
names(df.bed) <- c("chr_num", "start", "end", "peak", "maxp", "maxp_pos")
df.bed <- df.bed %>% mutate(chr = paste0("chr", chr_num)) %>% arrange(peak)

# Peaks 1/4/5 (chr6/16/17) = Bambomo (BAM) origin; 2/3/6 (chr8/13/24) = other.
peak_meta <- tribble(
  ~peak, ~chr,    ~origin,
  1,     "chr6",  "BAM",
  2,     "chr8",  "other",
  3,     "chr13", "other",
  4,     "chr16", "BAM",
  5,     "chr17", "BAM",
  6,     "chr24", "other"
)
df.bed <- df.bed %>% left_join(peak_meta, by = c("peak", "chr"))

# ---------------------------------------------------------------------------
# 2. Per-peak SNP membership (chr, pos, peak). allele1 = GWAS BP-associated
#    allele, but per-SNP BP orientation is re-derived empirically below.
# ---------------------------------------------------------------------------
snps_path <- file.path(root, "output_data", "06_Association",
                       "PKINGS_ALL_WOB_EXCLUDED",
                       "PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS.txt")
df.snps <- read.table(snps_path, header = TRUE, sep = "\t") %>%
  transmute(chr, POS = as.integer(ps), allele1, allele0, log_p, peak = as.integer(peak))

# ---------------------------------------------------------------------------
# 3. Sample metadata (PN matches VCF sample IDs)
# ---------------------------------------------------------------------------
fish_path <- file.path(root, "input_data", "01_Terra", "data_model",
                       "fish_data_2024a.txt")
df.fish <- read.table(fish_path, header = TRUE, sep = "\t", quote = "") %>%
  transmute(PN = trimws(PN), POP = trimws(POP), Phenotype = trimws(Phenotype))

# Reference groups for defining canonical BP / TP haplotypes
bam_bp_ids   <- df.fish %>% filter(POP == "BAM", Phenotype == "BP") %>% pull(PN)
other_bp_ids <- df.fish %>% filter(POP %in% c("BIK", "BAVA", "BIR"),
                                   Phenotype == "BP") %>% pull(PN)
# Peak 6 is a BP2-only locus (Biroundou carries no BP-class copy), so its
# BP reference is Bikagala + Bavavela BP alone.
bp2_bp_ids   <- df.fish %>% filter(POP %in% c("BIK", "BAVA"), Phenotype == "BP") %>% pull(PN)
tp_ref_ids   <- df.fish %>% filter(POP %in% c("APA", "BEN", "DOV", "DOG"),
                                   Phenotype == "TP") %>% pull(PN)
bp_ref_ids <- function(pk, origin) {
  if (origin == "BAM") bam_bp_ids else if (pk == 6) bp2_bp_ids else other_bp_ids
}

cat(sprintf("Reference haplotype pools: BAM BP=%d, BIK/BAVA/BIR BP (Peaks 2, 3)=%d, BIK/BAVA BP (Peak 6)=%d, TP (APA/BEN/DOV/DOG)=%d\n",
            length(bam_bp_ids), length(other_bp_ids), length(bp2_bp_ids), length(tp_ref_ids)))

# ---------------------------------------------------------------------------
# 4. Load phased GT matrix (bcftools query -H output)
#    header: #[1]CHROM  [2]POS  [3]REF  [4]ALT  [5]SAMPLE:GT ...
# ---------------------------------------------------------------------------
gt_path <- file.path(outdir, "PKINGS_WOB_peak_snps.gt_phased.txt")
hdr <- strsplit(readLines(gt_path, n = 1), "\t")[[1]]
hdr <- sub("^#?\\[[0-9]+\\]", "", hdr)   # strip leading #[N]
hdr <- sub(":GT$", "", hdr)              # strip :GT suffix from sample cols
gt <- fread(gt_path, header = FALSE, skip = 1, sep = "\t",
            colClasses = "character")
setnames(gt, hdr)
gt[, POS := as.integer(trimws(POS))]
gt[, CHROM := trimws(CHROM)]

sample_cols <- setdiff(names(gt), c("CHROM", "POS", "REF", "ALT"))
cat(sprintf("Phased GT matrix: %d SNPs x %d samples\n", nrow(gt), length(sample_cols)))

# Split each phased GT "a|b" into two allele matrices (0/1 = REF/ALT dosage).
# Missing alleles ('.') -> NA.
parse_allele <- function(v, which_hap) {
  a <- tstrsplit(v, "\\|", fixed = FALSE)[[which_hap]]
  a[a == "."] <- NA
  as.integer(a)
}
hapA <- as.matrix(gt[, lapply(.SD, parse_allele, which_hap = 1), .SDcols = sample_cols])
hapB <- as.matrix(gt[, lapply(.SD, parse_allele, which_hap = 2), .SDcols = sample_cols])
rownames(hapA) <- rownames(hapB) <- paste(gt$CHROM, gt$POS, sep = "_")
# ALT dosage genotype (0/1/2) for index-SNP concordance + orientation
dose <- hapA + hapB

# Map each matrix SNP row to its peak
snp_key <- gt[, .(CHROM, POS)][df.snps, on = c("CHROM" = "chr", "POS" = "POS"),
                               nomatch = 0]
# ensure alignment: build a lookup of row index -> peak
row_peak <- df.snps$peak[match(paste(gt$CHROM, gt$POS),
                               paste(df.snps$chr, df.snps$POS))]
stopifnot(!any(is.na(row_peak)))   # every matrix SNP must belong to a peak

# ---------------------------------------------------------------------------
# 5. Per-peak haplotype-allele calling
# ---------------------------------------------------------------------------
call_peak <- function(pk) {
  chr_pk <- peak_meta$chr[peak_meta$peak == pk]
  origin <- peak_meta$origin[peak_meta$peak == pk]
  bp_ref <- intersect(bp_ref_ids(pk, origin), sample_cols)
  tp_ref <- intersect(tp_ref_ids, sample_cols)

  rows <- which(row_peak == pk)
  snp_ids <- rownames(hapA)[rows]

  # Consensus (per-SNP major allele) among reference haplotypes.
  # Reference haplotypes = both copies of each reference individual.
  bp_hapmat <- cbind(hapA[rows, bp_ref, drop = FALSE], hapB[rows, bp_ref, drop = FALSE])
  tp_hapmat <- cbind(hapA[rows, tp_ref, drop = FALSE], hapB[rows, tp_ref, drop = FALSE])
  bp_cons <- as.integer(rowMeans(bp_hapmat, na.rm = TRUE) >= 0.5)
  tp_cons <- as.integer(rowMeans(tp_hapmat, na.rm = TRUE) >= 0.5)

  # Informative SNPs: consensus differs between BP and TP references.
  informative <- which(bp_cons != tp_cons)
  n_inf <- length(informative)

  # Classify a single haplotype vector against the two consensuses.
  classify_hap <- function(h) {
    if (n_inf == 0) return(NA_character_)
    hi <- h[informative]
    d_bp <- sum(hi != bp_cons[informative], na.rm = TRUE)
    d_tp <- sum(hi != tp_cons[informative], na.rm = TRUE)
    if (abs(d_bp - d_tp) <= AMBIG_MARGIN) return("ambig")
    if (d_bp < d_tp) "BP" else "TP"
  }

  hA_class <- vapply(sample_cols, function(s) classify_hap(hapA[rows, s]), character(1))
  hB_class <- vapply(sample_cols, function(s) classify_hap(hapB[rows, s]), character(1))

  bp_hap_count <- (hA_class == "BP") + (hB_class == "BP")
  peak_state <- dplyr::case_when(
    hA_class == "ambig" | hB_class == "ambig" |
      is.na(hA_class) | is.na(hB_class)          ~ "ambiguous",
    bp_hap_count == 2                            ~ "BP_hom",
    bp_hap_count == 1                            ~ "het",
    bp_hap_count == 0                            ~ "TP_hom"
  )

  states <- tibble(
    PN = sample_cols, peak = pk, chr = chr_pk, origin = origin,
    hapA_class = hA_class, hapB_class = hB_class,
    bp_hap_count = bp_hap_count, peak_state = peak_state,
    n_snps = length(rows), n_informative = n_inf
  )

  key <- tibble(
    peak = pk, chr = chr_pk, origin = origin,
    n_snps = length(rows), n_informative = n_inf,
    n_bp_ref = length(bp_ref), n_tp_ref = length(tp_ref),
    bp_consensus = paste0(bp_cons, collapse = ""),
    tp_consensus = paste0(tp_cons, collapse = "")
  )
  list(states = states, key = key)
}

res <- lapply(sort(unique(row_peak)), call_peak)
states <- bind_rows(lapply(res, `[[`, "states"))
allele_key <- bind_rows(lapply(res, `[[`, "key"))

# Attach phenotype metadata
states <- states %>% left_join(df.fish, by = "PN")

# ---------------------------------------------------------------------------
# 6. Index-SNP dosage concordance (robustness cross-check)
#    Reproduce the index-SNP dosage state from the same phased matrix and
#    compare to the multi-SNP haplotype state.
# ---------------------------------------------------------------------------
index_rows <- df.snps %>%
  group_by(peak) %>% slice_max(log_p, n = 1, with_ties = FALSE) %>% ungroup()

idx_state <- purrr::map_dfr(seq_len(nrow(index_rows)), function(i) {
  pk  <- index_rows$peak[i]
  rid <- paste(index_rows$chr[i], index_rows$POS[i], sep = "_")
  origin <- peak_meta$origin[peak_meta$peak == pk]
  bp_ref <- intersect(bp_ref_ids(pk, origin), sample_cols)
  tp_ref <- intersect(tp_ref_ids, sample_cols)
  d <- dose[rid, ]
  # BP allele = ALT if ALT more common in BP ref than TP ref, else REF
  bp_is_alt <- mean(d[bp_ref], na.rm = TRUE) > mean(d[tp_ref], na.rm = TRUE)
  bp_dosage <- if (bp_is_alt) d else 2L - d
  tibble(PN = sample_cols, peak = pk,
         index_state = dplyr::case_when(
           bp_dosage == 2 ~ "BP_hom", bp_dosage == 1 ~ "het",
           bp_dosage == 0 ~ "TP_hom", TRUE ~ "ambiguous"))
})

concord <- states %>%
  select(PN, peak, peak_state) %>%
  left_join(idx_state, by = c("PN", "peak")) %>%
  mutate(agree = peak_state == index_state)

concord_summary <- concord %>%
  group_by(peak) %>%
  summarise(n = n(),
            n_agree = sum(agree, na.rm = TRUE),
            pct_agree = round(100 * mean(agree, na.rm = TRUE), 1),
            .groups = "drop")
overall_agree <- round(100 * mean(concord$agree, na.rm = TRUE), 1)

# ---------------------------------------------------------------------------
# 7. Write outputs + sanity report
# ---------------------------------------------------------------------------
write.csv(states, file.path(outdir, "peak_haplotype_states.csv"), row.names = FALSE)
write.csv(allele_key, file.path(outdir, "peak_haplotype_allele_key.csv"), row.names = FALSE)
write.csv(concord_summary, file.path(outdir, "peak_haplotype_concordance.csv"), row.names = FALSE)

cat("\n=== Peak haplotype-allele key ===\n")
print(as.data.frame(allele_key %>% select(peak, chr, origin, n_snps, n_informative,
                                          n_bp_ref, n_tp_ref)))

cat("\n=== Sanity: mean BP-hap count by reference group and peak ===\n")
sanity <- states %>%
  mutate(refgrp = case_when(
    PN %in% bam_bp_ids   ~ "BAM_BP",
    PN %in% other_bp_ids ~ "otherBP",
    PN %in% tp_ref_ids   ~ "TP_ref",
    TRUE ~ "focal")) %>%
  filter(refgrp != "focal") %>%
  group_by(peak, origin, refgrp) %>%
  summarise(mean_bp_hap = round(mean(bp_hap_count, na.rm = TRUE), 2), .groups = "drop") %>%
  pivot_wider(names_from = refgrp, values_from = mean_bp_hap)
print(as.data.frame(sanity))

cat(sprintf("\n=== Index-SNP dosage concordance: overall %.1f%% ===\n", overall_agree))
print(as.data.frame(concord_summary))

cat(sprintf("\nWrote:\n  %s\n  %s\n  %s\n",
            file.path(outdir, "peak_haplotype_states.csv"),
            file.path(outdir, "peak_haplotype_allele_key.csv"),
            file.path(outdir, "peak_haplotype_concordance.csv")))
