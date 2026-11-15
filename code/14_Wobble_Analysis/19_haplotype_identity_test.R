#!/usr/bin/env Rscript
# ============================================================================
# 19_haplotype_identity_test.R
#
# Is a BP-class haplotype found outside its home populations the SAME haplotype
# as one in the home populations, or a similar one that can be told apart?
# (Home: Bambomo BP for Peaks 1, 4, 5; Bikagala/Bavavela/Biroundou BP for
# Peaks 2, 3; Bikagala/Bavavela BP for Peak 6, a BP2-only locus; for
# Biroundou at Peaks 2, 3, 6 the reference is Bikagala/Bavavela.)
#
# For every BP-class copy outside the home populations (from
# 17_peak_hap_match_audit.R, Fig. S4 SNP set), over all window SNPs:
#   nn_mismatch     mismatches to the nearest home BP-class copy
#   n_nearest       how many home copies are at that distance
#   novel_alleles   SNPs where the copy carries an allele that no home
#                   BP-class copy carries (a difference noise within the home
#                   pool cannot explain)
#   novel_detail    position, raw GT, AD, GQ of the fish at each novel allele
#   diff_detail     the same for each SNP where it differs from its nearest
#                   home copy, with how many home copies carry its allele
# Calibration per peak: for every home copy, mismatches to its nearest other
# home copy, and how many home copies have no identical twin.
#
# Inputs:
#   PKINGS_WOB_MULTILOCUS/peak_hap_zoom/figs4_haplotypes.rds (from 17)
#   output_data/08_Peak_Analysis/PKINGS_ALL_WOB_EXCLUDED/..._SNPS_IN_PEAKS_TOP250.vcf.gz
#     (Fig. S4 SNP set; raw GT:AD:DP:GQ)
#   output_data/07_Phasing/...phased.renamed.vcf.gz
# Outputs (peak_hap_zoom/):
#   identity_test_figS4snps.tsv, identity_calibration_figS4snps.tsv
#
# Run:  Rscript code/14_Wobble_Analysis/19_haplotype_identity_test.R   (needs bcftools)
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

root <- rprojroot::find_root(rprojroot::has_file(".git") | rprojroot::is_git_root)
setwd(root)
outdir <- file.path(root, "output_data", "14_Wobble_Analysis", "PKINGS_WOB_MULTILOCUS", "peak_hap_zoom")
dat <- readRDS(file.path(outdir, "figs4_haplotypes.rds"))
figs4 <- dat$figs4; bed <- dat$bed

figs4_vcf <- file.path(root, "output_data", "08_Peak_Analysis", "PKINGS_ALL_WOB_EXCLUDED",
                       "PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS_TOP250.vcf.gz")
phased_vcf <- file.path(root, "output_data", "07_Phasing",
                        "MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_SNPS_ONLY.phased.renamed.vcf.gz")

query <- function(cmd) {
  f <- tempfile(fileext = ".tsv"); system(paste(cmd, ">", f))
  h <- sub(":GT.*$", "", sub("^#?\\s*\\[[0-9]+\\]", "", strsplit(readLines(f, n = 1), "\t")[[1]]))
  d <- fread(f, header = FALSE, skip = 1, sep = "\t", colClasses = "character"); setnames(d, h)
  d[, POS := as.integer(POS)]; d
}
sites_f <- tempfile(fileext = ".tsv")
system(sprintf("bcftools query -f '%%CHROM\\t%%POS\\n' %s > %s", figs4_vcf, sites_f))
ph <- query(sprintf("bcftools view -T %s -R %s %s 2>/dev/null | bcftools query -H -f '%%CHROM\\t%%POS[\\t%%GT]\\n'",
                    sites_f, sites_f, phased_vcf))
rw <- query(sprintf("bcftools query -H -f '%%CHROM\\t%%POS[\\t%%GT:%%AD:%%DP:%%GQ]\\n' %s", figs4_vcf))
rw <- rw[match(paste(ph$CHROM, ph$POS), paste(rw$CHROM, rw$POS))]
stopifnot(!anyNA(rw$POS))

samples <- intersect(setdiff(names(ph), c("CHROM", "POS")), unique(figs4$PN))
allele <- function(v, w) { a <- tstrsplit(v, "[|/]")[[w]]; a[a == "."] <- NA; as.integer(a) }
H <- cbind(as.matrix(ph[, lapply(.SD, allele, w = 1), .SDcols = samples]),
           as.matrix(ph[, lapply(.SD, allele, w = 2), .SDcols = samples]))
colnames(H) <- c(paste0(samples, "|A"), paste0(samples, "|B"))
row_peak <- bed$peak[match(ph$CHROM, paste0("chr", bed$chr_num))]

raw_at <- function(i, pn) {
  x <- strsplit(rw[[pn]][i], ":")[[1]]
  sprintf("%d[%s AD=%s GQ=%s]", ph$POS[i], x[1], x[2], x[4])
}

home_pool <- function(pk, POP) {
  origin <- unique(figs4$origin[figs4$peak == pk])
  pops <- if (origin == "BAM") "BAM" else if (POP == "BIR" || pk == 6) c("BIK", "BAVA") else c("BIK", "BAVA", "BIR")
  figs4 %>% filter(peak == pk, POP %in% pops, Phenotype == "BP", hap_class == "BP") %>% pull(id)
}

# Candidates: BP-class copies outside the home populations, plus Biroundou
# copies at Peaks 2, 3, 6 (tested against Bikagala/Bavavela).
cands <- figs4 %>%
  filter(hap_class == "BP",
         (origin == "BAM" & POP != "BAM") |
         (origin == "other" & !(POP %in% c("BIK", "BAVA"))))

res <- bind_rows(lapply(seq_len(nrow(cands)), function(j) {
  cd <- cands[j, ]; pk <- cd$peak; h <- cd$id
  r <- which(row_peak == pk); home <- setdiff(home_pool(pk, cd$POP), h)
  Hp <- H[r, , drop = FALSE]
  d <- colSums(Hp[, home, drop = FALSE] != Hp[, h], na.rm = TRUE)
  nn <- names(which.min(d))
  novel <- which(vapply(seq_along(r), function(i) !(Hp[i, h] %in% Hp[i, home]), logical(1)))
  diffs <- which(Hp[, h] != Hp[, nn])
  tibble(peak = pk, id = h, PN = cd$PN, POP = cd$POP, Phenotype = cd$Phenotype,
         reference = if (cd$origin == "BAM") "BAM" else if (cd$POP == "BIR" || pk == 6) "BIK+BAVA" else "BIK+BAVA+BIR",
         n_snps = length(r), nn_mismatch = min(d), nearest = nn, n_nearest = sum(d == min(d)),
         novel_alleles = length(novel),
         novel_detail = paste(vapply(r[novel], raw_at, "", pn = cd$PN), collapse = "; "),
         diff_detail = paste(vapply(diffs, function(i) sprintf("%s on %d/%d home copies", raw_at(r[i], cd$PN),
                                                               sum(Hp[i, home] == Hp[i, h]), length(home)), ""),
                             collapse = "; "))
}))
write_tsv(res, file.path(outdir, "identity_test_figS4snps.tsv"))

cal <- bind_rows(lapply(sort(unique(row_peak)), function(pk) {
  r <- which(row_peak == pk)
  for_pop <- if (unique(figs4$origin[figs4$peak == pk]) == "BAM") "BAM" else "BIK"
  home <- home_pool(pk, for_pop)
  nn <- vapply(home, function(h) min(colSums(H[r, setdiff(home, h), drop = FALSE] != H[r, h])), numeric(1))
  tibble(peak = pk, reference = if (for_pop == "BAM") "BAM" else if (pk == 6) "BIK+BAVA" else "BIK+BAVA+BIR", n_snps = length(r),
         home_copies = length(home), nn_median = median(nn), nn_q95 = quantile(nn, 0.95), nn_max = max(nn),
         copies_without_identical_twin = sum(nn > 0))
}))
write_tsv(cal, file.path(outdir, "identity_calibration_figS4snps.tsv"))

cat("=== Calibration: home BP-class copies vs their nearest other home copy ===\n"); print(as.data.frame(cal))
cat("\n=== Candidates (P. kingsleyae; other species summarised by min mismatch) ===\n")
print(res %>% filter(!startsWith(POP, "P")) %>%
        select(peak, PN, POP, Phenotype, reference, n_snps, nn_mismatch, n_nearest, novel_alleles) %>%
        arrange(peak, nn_mismatch) %>% as.data.frame())
print(res %>% filter(startsWith(POP, "P")) %>% group_by(peak, POP) %>%
        summarise(copies = n(), min_mismatch = min(nn_mismatch), min_novel = min(novel_alleles), .groups = "drop") %>%
        as.data.frame())
