#!/usr/bin/env Rscript
# ============================================================================
# 11_power_sim_make.R
#
# Power simulation for the Bambomo wobble scans: what single-locus effect would
# the flat BAM GWAS have detected? Plants one causal SNP (real BAM genotypes,
# drawn from the LD-pruned GEMMA geno) and simulates a phenotype from it.
# power_sim.sb then tests each planted SNP with GEMMA (Wald, real kinship), and
# 11_power_sim_summarize.R scores p against the genome-wide permutation
# thresholds of the real scans (08_run_p0_gemma.sh).
#
# Designs (fish sets match BAM_P0.pheno columns, 07_make_p0_phenos.R):
#   CLS        16 cases / 42 controls (CLS_BIN fish). Case set drawn by weighted
#              sampling without replacement, weight = OR^dosage.
#   MANUAL     the 66 MANUAL_BIN fish: 16 true cases drawn as above, then 9
#              non-causal "manual wobbles" added at random -> 25 cases, the
#              observed manual/classifier disagreement (label noise).
#   MANUAL_CLEAN  66 fish, 25 true cases, no noise (upper bound for manual n).
#   QUANT      the 66 VP0_INT fish: y = sqrt(PVE) * z(dosage) + sqrt(1 - PVE) * e.
# One SNP per replicate, reused across effect sizes (paired), MAF binned within
# each design's fish set.
#
# Inputs:  arg1 = BAM_ONLY_BP0_TPWOB1.prune.for_gemma.geno
#          (output_data/06_Association/BAM_ONLY_BP0_TPWOB1/, from 06_Association/02_run_rel_mat.sh)
#          input_data/14_Wobble_Analysis/BAM_P0/BAM_P0.samples.tsv
# Outputs (output_data/14_Wobble_Analysis/BAM_P0/POWER_SIM/):
#   sim_design.tsv   sim_id, design, maf_bin, effect, rep, snp, maf
#   sim_geno.txt     one BIMBAM line per sim (rs = sim_id), 67 fish in GEMMA order
#   sim_pheno.tsv    67 rows x n_sim columns, "NA" = not in the design
#
# Run:  Rscript code/14_Wobble_Analysis/11_power_sim_make.R <prune.for_gemma.geno>
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

root <- rprojroot::find_root(rprojroot::is_git_root)
setwd(root)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1, file.exists(args[1]))

outdir <- file.path("output_data", "14_Wobble_Analysis", "BAM_P0", "POWER_SIM")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

set.seed(20261012)
N_REP     <- 100
ORS       <- c(1, 1.5, 2, 3, 5, 8, 12, 20, 50)   # 1 = null calibration, 50 = positive control
PVES      <- c(0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40)
MAF_BINS  <- c(0.05, 0.15, 0.30, 0.50)
N_TRUE    <- 16   # classifier-confident wobbles (CLS_BIN cases)
N_MANUAL  <- 25   # manual wobbles

samples <- read_tsv("input_data/14_Wobble_Analysis/BAM_P0/BAM_P0.samples.tsv", show_col_types = FALSE)
n_fish <- nrow(samples)

geno <- fread(args[1], header = FALSE)
stopifnot(ncol(geno) == 3 + n_fish)
G <- as.matrix(geno[, -(1:3)])
snp_id <- geno[[1]]

sets <- list(
  CLS          = which(!is.na(samples$CLS_BIN)),
  MANUAL       = which(!is.na(samples$MANUAL_BIN)),
  MANUAL_CLEAN = which(!is.na(samples$MANUAL_BIN)),
  QUANT        = which(!is.na(samples$VP0_INT))
)
stopifnot(length(sets$CLS) == 58, length(sets$MANUAL) == 66, length(sets$QUANT) == 66)

maf_in <- function(idx) {
  af <- rowMeans(G[, idx, drop = FALSE]) / 2
  pmin(af, 1 - af)
}

draw_cases <- function(g, n_case, or) {
  sample(seq_along(g), n_case, prob = or^g)
}

sim_rows <- list(); pheno_cols <- list(); geno_rows <- list(); k <- 0

for (design in names(sets)) {
  idx <- sets[[design]]
  maf <- maf_in(idx)
  effects <- if (design == "QUANT") PVES else ORS
  for (b in seq_len(length(MAF_BINS) - 1)) {
    pool <- which(maf >= MAF_BINS[b] & maf < MAF_BINS[b + 1] + (b == length(MAF_BINS) - 1) * 1e-9)
    snps <- sample(pool, N_REP)
    for (r in seq_len(N_REP)) {
      s <- snps[r]
      g <- G[s, idx]
      for (eff in effects) {
        y <- rep(NA_real_, n_fish)
        if (design == "QUANT") {
          z <- as.numeric(scale(g))
          y[idx] <- sqrt(eff) * z + sqrt(1 - eff) * rnorm(length(idx))
        } else {
          n_case <- if (design == "CLS") N_TRUE else if (design == "MANUAL") N_TRUE else N_MANUAL
          case <- draw_cases(g, n_case, eff)
          if (design == "MANUAL") {
            case <- c(case, sample(setdiff(seq_along(idx), case), N_MANUAL - N_TRUE))
          }
          yy <- rep(0, length(idx)); yy[case] <- 1
          y[idx] <- yy
        }
        k <- k + 1
        sid <- sprintf("sim%05d", k)
        sim_rows[[k]] <- tibble(sim_id = sid, design = design,
                                maf_bin = sprintf("%.2f-%.2f", MAF_BINS[b], MAF_BINS[b + 1]),
                                effect = eff, rep = r, snp = snp_id[s], maf = maf[s])
        pheno_cols[[k]] <- y
        geno_rows[[k]] <- paste(c(paste0(sid, ":", snp_id[s]), geno[[2]][s], geno[[3]][s], G[s, ]),
                                collapse = " ")
      }
    }
  }
}

design_tbl <- bind_rows(sim_rows)
write_tsv(design_tbl, file.path(outdir, "sim_design.tsv"))
writeLines(unlist(geno_rows), file.path(outdir, "sim_geno.txt"))
P <- do.call(cbind, pheno_cols)
P_chr <- ifelse(is.na(P), "NA", formatC(P, digits = 6, format = "g"))
fwrite(as.data.table(P_chr), file.path(outdir, "sim_pheno.tsv"), sep = "\t", col.names = FALSE, quote = FALSE)

design_tbl %>% count(design, maf_bin) %>% print(n = Inf)
cat("total sims:", nrow(design_tbl), "\n")
