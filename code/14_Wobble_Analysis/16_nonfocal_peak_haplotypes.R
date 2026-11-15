#!/usr/bin/env Rscript
# ============================================================================
# 16_nonfocal_peak_haplotypes.R
#
# Do fish OUTSIDE the focal BP/TP populations carry the BP-class HAPLOTYPE at
# each of the six GWAS peaks, or only the BP allele at the index SNP?
#
# Scores every phased haplotype of every fish in the phased VCF (305 fish: the
# 214 focal fish plus MOV, IVI, COB and the other Paramormyrops species)
# against the SAME per-peak BP and TP consensus haplotypes and diagnostic SNPs
# that 01_call_peak_haplotypes.R defined from the focal reference fish. The
# consensus strings are read from peak_haplotype_allele_key.csv; the script
# stops if the focal fish in the new extraction do not reproduce them, or do
# not reproduce the per-haplotype classes in peak_haplotype_states.csv.
#
# Per haplotype: fraction of diagnostic SNPs matching the BP consensus
# (missing excluded), class under the 01 rule (closer consensus; tie = ambig),
# BP allele at the index SNP (01 orientation rule), and the longest run of
# consecutive BP-matching diagnostic SNPs. Per fish (phase-free checks): the
# fraction of diagnostic SNPs homozygous for the BP allele, and, if an unphased
# VCF is present, the fraction of diagnostic SNPs with a raw call, phased-vs-raw
# genotype concordance and median depth (how much of the phased genotype was
# imputed by SHAPEIT4).
#
# Inputs (all relative to repo root):
#   phased VCF: $phased_VCF (output_data/07_Phasing/...phased.renamed.vcf.gz)
#   unphased VCF (optional QC): output_data/03_QC/...FILTERED_BIALLELIC_SNPS_ONLY.renamed.vcf.gz or
#     ...FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz (03_VCF_QC_Filter)
#   output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS.txt
#   output_data/14_Wobble_Analysis/PKINGS_WOB_MULTILOCUS/peak_haplotype_allele_key.csv
#   output_data/14_Wobble_Analysis/PKINGS_WOB_MULTILOCUS/peak_haplotype_states.csv     (check)
#   output_data/14_Wobble_Analysis/PKINGS_WOB_MULTILOCUS/PKINGS_WOB_peak_snps.gt_phased.txt (check)
#   input_data/01_Terra/data_model/fish_data_2024a.txt
#
# Outputs (PKINGS_WOB_MULTILOCUS/nonfocal_haplotypes/):
#   per_haplotype.tsv           one row per fish x peak x haplotype copy
#   summary_by_group_peak.tsv   group x peak: n, median/range BP-match, n >= 0.90, n BP class
#
# Needs bcftools on PATH. Runs locally in seconds (215 SNPs).
# Run:  Rscript code/14_Wobble_Analysis/16_nonfocal_peak_haplotypes.R
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

root <- rprojroot::find_root(rprojroot::has_file(".git") | rprojroot::is_git_root)
setwd(root)

mldir  <- file.path("output_data", "14_Wobble_Analysis", "PKINGS_WOB_MULTILOCUS")
outdir <- file.path(mldir, "nonfocal_haplotypes")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

NEAR_COMPLETE <- 0.90   # BP-match fraction for a "near-complete BP haplotype"
AMBIG_MARGIN  <- 0      # as in 01_call_peak_haplotypes.R

vcf_name <- "MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_SNPS_ONLY.phased.renamed.vcf.gz"
phased_vcf <- c(file.path("output_data", "07_Phasing", vcf_name))
phased_vcf <- phased_vcf[file.exists(phased_vcf)][1]
stopifnot(!is.na(phased_vcf))
raw_vcf <- c(file.path("output_data", "03_QC",
                       "MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_SNPS_ONLY.renamed.vcf.gz"),
             file.path("output_data", "03_QC",
                       "MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz"))
raw_vcf <- raw_vcf[file.exists(raw_vcf)][1]
cat(sprintf("Phased VCF: %s\nUnphased VCF (QC): %s\n", phased_vcf,
            ifelse(is.na(raw_vcf), "none found, raw-call QC skipped", raw_vcf)))

# ---------------------------------------------------------------------------
# 1. Peak definitions, peak SNPs, metadata (as in 01_call_peak_haplotypes.R)
# ---------------------------------------------------------------------------
peak_meta <- tribble(
  ~peak, ~chr,    ~origin,
  1,     "chr6",  "BAM",
  2,     "chr8",  "other",
  3,     "chr13", "other",
  4,     "chr16", "BAM",
  5,     "chr17", "BAM",
  6,     "chr24", "other"
)

df.snps <- read.table(file.path("output_data", "06_Association", "PKINGS_ALL_WOB_EXCLUDED",
                                "PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS.txt"),
                      header = TRUE, sep = "\t") %>%
  transmute(chr, POS = as.integer(ps), log_p, peak = as.integer(peak))

df.fish <- read.table(file.path("input_data", "01_Terra", "data_model", "fish_data_2024a.txt"),
                      header = TRUE, sep = "\t", quote = "") %>%
  transmute(PN = trimws(PN), POP = trimws(POP), Phenotype = trimws(Phenotype))

bam_bp_ids   <- df.fish %>% filter(POP == "BAM", Phenotype == "BP") %>% pull(PN)
other_bp_ids <- df.fish %>% filter(POP %in% c("BIK", "BAVA", "BIR"), Phenotype == "BP") %>% pull(PN)
# Peak 6 is a BP2-only locus: its BP reference is Bikagala + Bavavela alone, as in 01.
bp2_bp_ids   <- df.fish %>% filter(POP %in% c("BIK", "BAVA"), Phenotype == "BP") %>% pull(PN)
tp_ref_ids   <- df.fish %>% filter(POP %in% c("APA", "BEN", "DOV", "DOG"), Phenotype == "TP") %>% pull(PN)
focal_pops   <- c("BAM", "APA", "BEN", "BIK", "BAVA", "BIR", "DOV", "DOG")

# ---------------------------------------------------------------------------
# 2. Extract phased (and optionally raw) genotypes at the peak SNPs, all fish
# ---------------------------------------------------------------------------
regions <- tempfile(fileext = ".txt")
write.table(df.snps[, c("chr", "POS")], regions, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

read_query <- function(vcf, fmt) {
  tmp <- tempfile(fileext = ".txt")
  status <- system2("bcftools", c("query", "-H", "-R", regions, "-f", shQuote(fmt), vcf),
                    stdout = tmp, stderr = FALSE)
  stopifnot(status == 0)
  hdr <- strsplit(readLines(tmp, n = 1), "\t")[[1]]
  hdr <- sub("^#?\\[[0-9]+\\]", "", hdr)
  hdr <- sub(":.*$", "", hdr)
  d <- fread(tmp, header = FALSE, skip = 1, sep = "\t", colClasses = "character")
  setnames(d, hdr)
  d[, POS := as.integer(POS)]
  unlink(tmp)
  d
}

gt <- read_query(phased_vcf, "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n")
sample_cols <- setdiff(names(gt), c("CHROM", "POS", "REF", "ALT"))
cat(sprintf("Phased GT matrix: %d SNPs x %d samples\n", nrow(gt), length(sample_cols)))

row_peak <- df.snps$peak[match(paste(gt$CHROM, gt$POS), paste(df.snps$chr, df.snps$POS))]
stopifnot(!any(is.na(row_peak)), nrow(gt) == nrow(df.snps))

parse_allele <- function(v, which_hap) {
  a <- tstrsplit(v, "|", fixed = TRUE)[[which_hap]]
  a[a == "."] <- NA
  as.integer(a)
}
hapA <- as.matrix(gt[, lapply(.SD, parse_allele, which_hap = 1), .SDcols = sample_cols])
hapB <- as.matrix(gt[, lapply(.SD, parse_allele, which_hap = 2), .SDcols = sample_cols])
dose <- hapA + hapB

# Check 1: focal fish genotypes identical to the matrix 01 used.
old <- fread(file.path(mldir, "PKINGS_WOB_peak_snps.gt_phased.txt"), header = FALSE, skip = 1,
             sep = "\t", colClasses = "character")
old_hdr <- sub(":GT$", "", sub("^#?\\[[0-9]+\\]", "",
                               strsplit(readLines(file.path(mldir, "PKINGS_WOB_peak_snps.gt_phased.txt"),
                                                  n = 1), "\t")[[1]]))
setnames(old, old_hdr)
old[, POS := as.integer(POS)]
stopifnot(identical(old$POS, gt$POS))
focal_in_old <- setdiff(names(old), c("CHROM", "POS", "REF", "ALT"))
stopifnot(all(focal_in_old %in% sample_cols))
stopifnot(all(as.matrix(old[, ..focal_in_old]) == as.matrix(gt[, ..focal_in_old])))
cat(sprintf("Check 1 passed: %d focal fish have identical phased genotypes to the 01 input.\n",
            length(focal_in_old)))

# ---------------------------------------------------------------------------
# 3. Consensus haplotypes + diagnostic SNPs, from the 01 allele key
# ---------------------------------------------------------------------------
allele_key <- read.csv(file.path(mldir, "peak_haplotype_allele_key.csv"),
                       colClasses = c(bp_consensus = "character", tp_consensus = "character"))

old_states <- read.csv(file.path(mldir, "peak_haplotype_states.csv"))

# Index SNP = top log_p SNP per peak (as in 01).
index_rows <- df.snps %>% group_by(peak) %>%
  slice_max(log_p, n = 1, with_ties = FALSE) %>% ungroup()

longest_run <- function(match, pos) {
  # longest run of consecutive diagnostic SNPs matching BP (NA skipped)
  ok <- !is.na(match); m <- match[ok]; p <- pos[ok]
  if (!any(m)) return(c(n = 0, bp = 0))
  r <- rle(m); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  w <- which(r$values)
  w <- w[which.max(r$lengths[w])]
  c(n = r$lengths[w], bp = p[ends[w]] - p[starts[w]])
}

score_peak <- function(pk) {
  key    <- allele_key[allele_key$peak == pk, ]
  origin <- key$origin
  rows   <- which(row_peak == pk)
  pos    <- gt$POS[rows]
  bp_cons <- as.integer(strsplit(key$bp_consensus, "")[[1]])
  tp_cons <- as.integer(strsplit(key$tp_consensus, "")[[1]])
  stopifnot(length(bp_cons) == length(rows))

  # Check 2: the focal reference fish reproduce the stored consensus.
  bp_ref <- intersect(if (origin == "BAM") bam_bp_ids else if (pk == 6) bp2_bp_ids else other_bp_ids, sample_cols)
  tp_ref <- intersect(tp_ref_ids, sample_cols)
  re_bp <- as.integer(rowMeans(cbind(hapA[rows, bp_ref], hapB[rows, bp_ref]), na.rm = TRUE) >= 0.5)
  re_tp <- as.integer(rowMeans(cbind(hapA[rows, tp_ref], hapB[rows, tp_ref]), na.rm = TRUE) >= 0.5)
  stopifnot(identical(re_bp, bp_cons), identical(re_tp, tp_cons))

  inf <- which(bp_cons != tp_cons)
  stopifnot(length(inf) == key$n_informative)
  pos_inf <- pos[inf]

  # BP allele at the index SNP (01 rule: ALT if ALT commoner in BP ref than TP ref)
  irow <- which(gt$CHROM == index_rows$chr[index_rows$peak == pk] &
                gt$POS == index_rows$POS[index_rows$peak == pk])
  bp_is_alt <- mean(dose[irow, bp_ref], na.rm = TRUE) > mean(dose[irow, tp_ref], na.rm = TRUE)
  bp_allele_idx <- as.integer(bp_is_alt)

  one_hap <- function(H, s, label) {
    h <- H[rows, s][inf]
    match <- h == bp_cons[inf]
    n_called <- sum(!is.na(match))
    n_bp <- sum(match, na.rm = TRUE)
    d_bp <- n_called - n_bp; d_tp <- n_bp           # biallelic: mismatch BP = match TP
    cls <- if (n_called == 0) NA_character_ else if (abs(d_bp - d_tp) <= AMBIG_MARGIN) "ambig" else
      if (d_bp < d_tp) "BP" else "TP"
    lr <- longest_run(match, pos_inf)
    tibble(PN = s, peak = pk, hap = label, n_diag = length(inf), n_called = n_called,
           n_bp_match = n_bp, bp_match_frac = n_bp / n_called, hap_class = cls,
           index_pos = gt$POS[irow], index_bp_allele = as.integer(H[irow, s] == bp_allele_idx),
           longest_bp_run_snps = unname(lr["n"]), longest_bp_run_bp = unname(lr["bp"]))
  }
  haps <- bind_rows(lapply(sample_cols, function(s)
    bind_rows(one_hap(hapA, s, "A"), one_hap(hapB, s, "B"))))

  # Phase-free per fish: fraction of diagnostic SNPs homozygous for the BP allele.
  hom_bp <- colMeans(hapA[rows[inf], , drop = FALSE] == bp_cons[inf] &
                     hapB[rows[inf], , drop = FALSE] == bp_cons[inf], na.rm = TRUE)
  haps$fish_bp_hom_frac <- hom_bp[haps$PN]
  haps$chr <- key$chr; haps$origin <- origin
  haps$diag_first <- min(pos_inf); haps$diag_last <- max(pos_inf)
  haps
}

per_hap <- bind_rows(lapply(sort(unique(row_peak)), score_peak))

# Check 3: focal fish reproduce the 01 per-haplotype classes.
chk <- per_hap %>% select(PN, peak, hap, hap_class) %>%
  pivot_wider(names_from = hap, values_from = hap_class, names_prefix = "new_") %>%
  inner_join(old_states %>% select(PN, peak, hapA_class, hapB_class), by = c("PN", "peak"))
stopifnot(nrow(chk) == nrow(old_states),
          all(chk$new_A == chk$hapA_class), all(chk$new_B == chk$hapB_class))
cat(sprintf("Check 3 passed: %d focal fish x peak rows reproduce the 01 haplotype classes.\n", nrow(chk)))

# ---------------------------------------------------------------------------
# 4. Optional raw-call QC (unphased VCF): call rate, concordance, depth
# ---------------------------------------------------------------------------
if (!is.na(raw_vcf)) {
  raw <- read_query(raw_vcf, "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP]\n")
  # header after stripping: CHROM POS REF ALT then sample, sample (GT, DP) pairs
  gt_cols <- seq(5, ncol(raw), by = 2); dp_cols <- gt_cols + 1
  raw_samples <- names(raw)[gt_cols]
  stopifnot(identical(raw$POS, gt$POS), identical(raw$REF, gt$REF), identical(raw$ALT, gt$ALT),
            all(sample_cols %in% raw_samples))
  rgt <- as.matrix(raw[, ..gt_cols]); colnames(rgt) <- raw_samples
  rdp <- as.matrix(raw[, ..dp_cols]); colnames(rdp) <- raw_samples
  rdose <- matrix(NA_integer_, nrow(rgt), ncol(rgt), dimnames = dimnames(rgt))
  for (j in seq_len(ncol(rgt))) {
    g <- rgt[, j]
    ok <- !grepl(".", g, fixed = TRUE)
    rdose[ok, j] <- vapply(strsplit(g[ok], "[/|]"), function(x) sum(as.integer(x)), integer(1))
  }
  suppressWarnings(storage.mode(rdp) <- "integer")
  raw_qc <- map_dfr(sort(unique(row_peak)), function(pk) {
    key <- allele_key[allele_key$peak == pk, ]
    rows <- which(row_peak == pk)
    inf <- rows[strsplit(key$bp_consensus, "")[[1]] != strsplit(key$tp_consensus, "")[[1]]]
    tibble(PN = sample_cols, peak = pk,
           raw_called_frac = colMeans(!is.na(rdose[inf, sample_cols, drop = FALSE])),
           raw_phased_concord = colSums(rdose[inf, sample_cols, drop = FALSE] ==
                                        dose[inf, sample_cols, drop = FALSE], na.rm = TRUE) /
             colSums(!is.na(rdose[inf, sample_cols, drop = FALSE])),
           median_dp = apply(rdp[inf, sample_cols, drop = FALSE], 2, median, na.rm = TRUE))
  })
  per_hap <- per_hap %>% left_join(raw_qc, by = c("PN", "peak"))
} else {
  per_hap <- per_hap %>% mutate(raw_called_frac = NA_real_, raw_phased_concord = NA_real_,
                                median_dp = NA_real_)
}

# ---------------------------------------------------------------------------
# 5. Groups, outputs
# ---------------------------------------------------------------------------
per_hap <- per_hap %>%
  left_join(df.fish, by = "PN") %>%
  mutate(POP = coalesce(POP, sub("_.*$", "", PN)),
         focal = POP %in% focal_pops,
         group = case_when(focal ~ paste(POP, Phenotype, sep = "_"),
                           POP == "COB" ~ paste(POP, Phenotype, sep = "_"),   # BP Mengono / TP Iboundji
                           TRUE ~ POP))

group_order <- c("BAM_BP", "BAM_Wobble", "BAM_TP", "BIR_BP", "BIK_BP", "BIK_Wobble",
                 "BAVA_BP", "BAVA_Wobble", "APA_TP", "APA_Wobble", "BEN_TP", "DOV_TP", "DOG_TP",
                 "MOV", "IVI", "COB_BP", "COB_TP", "PCUR", "PLON", "PNTO", "PSN2", "PVAD",
                 "PPAR", "PSZA", "PGAB", "PHOP", "PMAG")
stopifnot(all(per_hap$group %in% group_order))

per_hap <- per_hap %>%
  mutate(group = factor(group, levels = group_order)) %>%
  arrange(peak, group, PN, hap) %>%
  mutate(across(c(bp_match_frac, fish_bp_hom_frac, raw_called_frac, raw_phased_concord),
                ~ round(.x, 3))) %>%
  select(PN, group, POP, Phenotype, focal, peak, chr, origin, hap, n_diag, n_called,
         n_bp_match, bp_match_frac, hap_class, index_pos, index_bp_allele,
         longest_bp_run_snps, longest_bp_run_bp, fish_bp_hom_frac,
         raw_called_frac, raw_phased_concord, median_dp)
write_tsv(per_hap, file.path(outdir, "per_haplotype.tsv"))

diag_span <- per_hap %>% distinct(peak) %>%
  left_join(bind_rows(lapply(sort(unique(row_peak)), function(pk) {
    key <- allele_key[allele_key$peak == pk, ]
    pos <- gt$POS[row_peak == pk]
    inf <- strsplit(key$bp_consensus, "")[[1]] != strsplit(key$tp_consensus, "")[[1]]
    tibble(peak = pk, diag_start = min(pos[inf]), diag_end = max(pos[inf]),
           diag_span_bp = max(pos[inf]) - min(pos[inf]))
  })), by = "peak")

summary_tbl <- per_hap %>%
  group_by(group, focal, peak, chr, origin) %>%
  summarise(n_fish = n_distinct(PN), n_hap = n(), n_diag = first(n_diag),
            median_bp_match = median(bp_match_frac, na.rm = TRUE),
            min_bp_match = min(bp_match_frac, na.rm = TRUE),
            max_bp_match = max(bp_match_frac, na.rm = TRUE),
            n_hap_ge90 = sum(bp_match_frac >= NEAR_COMPLETE, na.rm = TRUE),
            n_hap_BP_class = sum(hap_class == "BP", na.rm = TRUE),
            n_hap_ambig = sum(hap_class == "ambig", na.rm = TRUE),
            n_hap_index_bp = sum(index_bp_allele, na.rm = TRUE),
            median_longest_run_bp = median(longest_bp_run_bp, na.rm = TRUE),
            median_raw_called = median(raw_called_frac, na.rm = TRUE),
            median_dp = median(median_dp, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(diag_span, by = "peak") %>%
  mutate(across(c(median_bp_match, min_bp_match, max_bp_match, median_raw_called), ~ round(.x, 3))) %>%
  arrange(peak, group)
write_tsv(summary_tbl, file.path(outdir, "summary_by_group_peak.tsv"))

cat("\n=== Diagnostic SNP span per peak ===\n")
print(as.data.frame(diag_span %>% left_join(allele_key %>% select(peak, chr, n_snps, n_informative),
                                            by = "peak")))
cat("\n=== BP-match by group x peak (n_hap_ge90 / n_hap_BP_class / n_hap_index_bp of n_hap) ===\n")
print(as.data.frame(summary_tbl %>% select(group, peak, n_hap, median_bp_match, min_bp_match,
                                           max_bp_match, n_hap_ge90, n_hap_BP_class, n_hap_index_bp)),
      row.names = FALSE)
cat(sprintf("\nWrote:\n  %s\n  %s\n", file.path(outdir, "per_haplotype.tsv"),
            file.path(outdir, "summary_by_group_peak.tsv")))
