#!/usr/bin/env Rscript
## 07_phaseB.R
## Phase B -- Peak 6 against the fixed-difference (drift-made cluster) null.
## METHODS.md part B.  Uses the IDENTICAL operator as sweep_detection/evidence/02_region_stats.R.
suppressPackageStartupMessages({ library(data.table); library(rprojroot) })
## accept .git as a directory (clone) or a file (linked worktree)
root <- find_root(has_file(".git") | has_dir(".git"))
source(file.path(root, "code/11_Drift_vs_Selection/sweep_detection/evidence/region_stats_lib.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_of <- function(f, d) { i <- match(f, args); if (is.na(i) || i == length(args)) d else args[i + 1] }
CLUST <- arg_of("--cluster_dir", Sys.getenv("CLUSTER_DIR",
                 file.path(root, "output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree")))
OUTD  <- arg_of("--out_dir", Sys.getenv("FCLR_OUT",
                 file.path(root, "output_data/11_Drift_vs_Selection/folded_clr")))
if (!startsWith(OUTD, "/")) OUTD <- file.path(root, OUTD)
dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
MIN_CLUSTERS <- 200L

dir_ev <- file.path(root, "output_data/11_Drift_vs_Selection/sweep_detection/evidence")
stats  <- fread(file.path(dir_ev, "region_stats.tsv.gz"))

## ---------------------------------------------------------------- B1 precondition
counts <- fread(file.path(CLUST, "fixed_diff_cluster_counts.tsv"))
stopifnot(all(c("pair_id","threshold","direction","n_snps","n_clusters") %in% names(counts)))
pairs_want <- data.table(
  pair    = c("A1w", "A2w", "A3w"),
  pop     = c("BAM", "BIK", "BAVA"),
  lineage = c("BP1", "BP2", "BP2"),
  peak    = c(1L, 6L, 6L))

pre <- rbindlist(lapply(seq_len(nrow(pairs_want)), function(i) {
  pr <- pairs_want$pair[i]
  bed <- file.path(CLUST, sprintf("fixed_diff_clusters.%s.bed", pr))
  n_bed <- if (file.exists(bed)) length(readLines(bed)) else NA_integer_
  bed9 <- file.path(CLUST, sprintf("fixed_diff_clusters_0.9.%s.bed", pr))
  n_bed9 <- if (file.exists(bed9)) length(readLines(bed9)) else NA_integer_
  c95 <- counts[pair_id == pr & abs(threshold - 0.95) < 1e-9 & direction == "BP_derived"]
  c90 <- counts[pair_id == pr & abs(threshold - 0.90) < 1e-9 & direction == "BP_derived"]
  data.table(pair = pr, population = pairs_want$pop[i], lineage = pairs_want$lineage[i],
             peak = pairs_want$peak[i],
             n_clusters_0.95_table = c95$n_clusters, n_clusters_0.95_bed = n_bed,
             n_snps_0.95_table = c95$n_snps,
             n_clusters_0.90_table = c90$n_clusters, n_clusters_0.90_bed = n_bed9,
             table_bed_agree = identical(as.integer(c95$n_clusters), as.integer(n_bed)),
             min_required = MIN_CLUSTERS,
             precondition = if (isTRUE(c95$n_clusters >= MIN_CLUSTERS)) "PASS" else "FAIL")
}))
fwrite(pre, file.path(OUTD, "phaseB_precondition.tsv"), sep = "\t")
cat("--- B1 precondition ---\n"); print(pre)
if (any(!pre$table_bed_agree)) stop("cluster count table disagrees with the BED line count")

usable <- pre[precondition == "PASS", pair]
if (!all(c("A2w", "A3w") %in% usable))
  stop("Peak 6 lineages do not both clear the precondition; Phase B not run")

## ---------------------------------------------------------------- B2 null loci
read_clusters <- function(pr) {
  f <- file.path(CLUST, sprintf("fixed_diff_clusters.%s.bed", pr))
  if (length(readLines(f)) == 0L) return(NULL)
  d <- fread(f, col.names = c("chrom", "start", "end", "n_snps", "centre_pos", "delta"))
  d[, .(region_id = paste0(pr, "_", chrom, "_", centre_pos), region_type = "fdcluster",
        pair = pr, chrom = chrom, centre = as.numeric(centre_pos),
        cluster_n_snps = n_snps, cluster_delta = delta)]
}
cl <- rbindlist(lapply(c("A2w", "A3w"), read_clusters))
cat(sprintf("\nNull loci: %d clusters (A2w %d, A3w %d) on %d chromosomes\n",
            nrow(cl), sum(cl$pair == "A2w"), sum(cl$pair == "A3w"), uniqueN(cl$chrom)))

## ---------------------------------------------------------------- B3 the same operator
## settings are asserted against the operator table written by 02_region_stats.R
ops_tbl <- fread(file.path(dir_ev, "region_stats_operators.tsv"))
spec <- list(
  rsb        = list(min_sites = 20L, direction = "high", operator = "window_mean",
                    load = function() load_rsb(root, pair_of("BP2"))),
  clr        = list(min_sites =  3L, direction = "high", operator = "region_max",
                    load = function() load_clr(root, "BP2")),
  clr_folded = list(min_sites =  3L, direction = "high", operator = "region_max",
                    load = function() load_clr(root, "BP2", suffix = "genome_wide_clr_folded.txt")),
  pi_ratio   = list(min_sites =  4L, direction = "low",  operator = "window_mean",
                    load = function() load_pi_ratio(root, "BP2", "TP2")))

vals <- list(); summ <- list()
for (st in names(spec)) {
  s <- spec[[st]]
  op <- make_operator(half_width = 100e3, window = 20e3, step = 5e3,
                      min_sites = s$min_sites, label = st)
  ref <- ops_tbl[statistic == st & lineage == "BP2"]
  if (nrow(ref)) {
    ref_op <- list(half_width = ref$half_width[1], window = ref$window[1],
                   step = ref$step[1], min_sites = ref$min_sites[1])
    assert_same_operator(op, ref_op, sprintf("Phase B and 02_region_stats.R (%s)", st))
    stopifnot(identical(ref$direction[1], s$direction), identical(ref$operator[1], s$operator))
  } else stop("no operator row for ", st, " BP2")

  dt <- s$load()
  cat(sprintf("%-11s %s | %d scored sites genome-wide\n", st, format_operator(op), nrow(dt)))
  r <- apply_operator(dt, cl[, .(region_id, chrom, centre)], op, s$direction, s$operator)
  r[, statistic := st]
  vals[[st]] <- merge(r, cl, by = "region_id")

  peak_val <- stats[statistic == st & lineage == "BP2" & region_id == "peak6", value]
  if (!length(peak_val)) stop("no peak6 BP2 value for ", st)
  nv <- r$value[!is.na(r$value)]
  pct <- if (s$direction == "high") 100 * mean(nv <= peak_val) else 100 * mean(nv >= peak_val)
  summ[[st]] <- data.table(
    statistic = st, primary = st %in% c("rsb", "clr", "pi_ratio"), direction = s$direction,
    peak6_value = peak_val, n_clusters_total = nrow(r), n_clusters_scored = length(nv),
    n_clusters_dropped_min_sites = nrow(r) - length(nv),
    null_median = median(nv), null_q95 = as.numeric(quantile(nv, if (s$direction == "high") 0.95 else 0.05)),
    peak6_pct_in_fdnull = pct,
    exceeds_95th = pct >= 95)
  ## per-lineage split
  for (pr in c("A2w", "A3w")) {
    nvp <- vals[[st]][pair == pr & !is.na(value), value]
    pctp <- if (s$direction == "high") 100 * mean(nvp <= peak_val) else 100 * mean(nvp >= peak_val)
    summ[[paste0(st, "_", pr)]] <- data.table(
      statistic = st, primary = st %in% c("rsb", "clr", "pi_ratio"), direction = s$direction,
      peak6_value = peak_val, n_clusters_total = sum(vals[[st]]$pair == pr),
      n_clusters_scored = length(nvp), n_clusters_dropped_min_sites = sum(vals[[st]]$pair == pr) - length(nvp),
      null_median = median(nvp), null_q95 = as.numeric(quantile(nvp, if (s$direction == "high") 0.95 else 0.05)),
      peak6_pct_in_fdnull = pctp, exceeds_95th = pctp >= 95, null_set = pr)
  }
}
S <- rbindlist(summ, fill = TRUE)
S[is.na(null_set), null_set := "A2w+A3w pooled"]
setcolorder(S, c("statistic", "null_set"))
setorder(S, statistic, null_set)
fwrite(S, file.path(OUTD, "phaseB_percentiles.tsv"), sep = "\t")
fwrite(rbindlist(vals, fill = TRUE), file.path(OUTD, "phaseB_cluster_values.tsv.gz"), sep = "\t")
cat("\n--- B3/B4 Phase B percentiles ---\n"); print(S)

## ---------------------------------------------------------------- B4 verdict
for (nsset in unique(S$null_set)) {
  for (clr_use in c("clr", "clr_folded")) {
    prim <- S[null_set == nsset & statistic %in% c("rsb", clr_use, "pi_ratio")]
    n_ok <- sum(prim$exceeds_95th)
    cat(sprintf("null %-15s CLR=%-10s : %d of 3 primaries >= 95th -> %s\n",
                nsset, clr_use, n_ok, if (n_ok >= 2) "SURVIVES Phase B" else "does NOT survive Phase B"))
  }
}
cat("\nPeak 1 / BP1 (Bambomo, pair A1w): precondition FAIL, recorded in phaseB_precondition.tsv;",
    "no substitute null run.\n")
