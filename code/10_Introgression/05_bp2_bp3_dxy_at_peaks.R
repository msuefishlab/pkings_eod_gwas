## 05_bp2_bp3_dxy_at_peaks.R
## JRG
## BP2 (Bikagala, Bavavela) and BP3 (Biroundou) share BP-associated haplotypes at Peaks 2 and 3. Was that haplotype
## exchanged recently (gene flow after the lineages split) or is it older?
##
## Raw dxy is confounded by local diversity: the peaks are divergent for every population pair. So the test statistic is
## relative node depth, RND = dxy(BP2, BP3) / mean(dxy(BP2, out), dxy(BP3, out)) with P. szaboii (PSZA) as the outgroup,
## which cancels local mutation-rate and diversity differences. Recent gene flow lowers RND; an old shared haplotype
## does not. Raw dxy and the TP2-BP3 comparison are reported alongside.
##
## Uses the pixy 5-kb dxy windows (output_data/05_PopGen/..._dxy.txt) and the extended peak windows (input_data/10_Introgression/peaks.bed).
## Lineage dxy pools the per-population-pair counts: sum(count_diffs) / sum(count_comparisons).
## Null: non-overlapping blocks of the same number of consecutive 5-kb windows, within chromosomes, excluding windows
## that overlap any extended peak.
##
## Run from the repo root: Rscript code/10_Introgression/05_bp2_bp3_dxy_at_peaks.R

library(data.table)

root <- system("git rev-parse --show-toplevel", intern = TRUE)
dxy_file <- file.path(root, "output_data/05_PopGen/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites_dxy.txt")
peaks_file <- file.path(root, "input_data/10_Introgression/peaks.bed")
out_file <- file.path(root, "output_data/10_Introgression/bp2_bp3_dxy_peaks.tsv")

groups <- list(
  "BP2-BP3" = c("BAVA_BIR", "BIK_BIR"),
  "TP2-BP3" = c("BIR_DOG", "BIR_DOV"),
  "BP2-TP2" = c("BAVA_DOG", "BAVA_DOV", "BIK_DOG", "BIK_DOV"),
  "BP2-OUT" = c("BAVA_PSZA", "BIK_PSZA"),
  "BP3-OUT" = c("BIR_PSZA"),
  "TP2-OUT" = c("DOG_PSZA", "DOV_PSZA")
)
pair_to_group <- data.table(pair = unlist(groups), group = rep(names(groups), lengths(groups)))

dxy <- fread(dxy_file, select = c("pop1", "pop2", "chromosome", "window_pos_1", "window_pos_2",
                                  "count_diffs", "count_comparisons"))
dxy[, pair := paste(pop1, pop2, sep = "_")]
dxy <- merge(dxy[pair %in% pair_to_group$pair], pair_to_group, by = "pair")
win <- dxy[, .(diffs = sum(count_diffs, na.rm = TRUE), comps = sum(count_comparisons, na.rm = TRUE)),
           by = .(group, chr = chromosome, start = window_pos_1, end = window_pos_2)]
# wide: one row per window, diffs/comps per group
wide <- dcast(win, chr + start + end ~ group, value.var = c("diffs", "comps"))
setorder(wide, chr, start)

peaks <- fread(peaks_file, col.names = c("chr", "start", "end", "peak", "max_logp", "index_bp"))
peaks[, chr := paste0("chr", chr)]
wide[, in_peak := FALSE]
for (i in seq_len(nrow(peaks))) {
  wide[chr == peaks$chr[i] & start <= peaks$end[i] & end >= peaks$start[i], in_peak := TRUE]
}

# Summaries for a set of windows: dxy per group, RND for BP2-BP3 and TP2-BP3
summarise_windows <- function(w) {
  d <- function(g) sum(w[[paste0("diffs_", g)]], na.rm = TRUE) / sum(w[[paste0("comps_", g)]], na.rm = TRUE)
  data.table(dxy_BP2_BP3 = d("BP2-BP3"), dxy_TP2_BP3 = d("TP2-BP3"), dxy_BP2_TP2 = d("BP2-TP2"),
             rnd_BP2_BP3 = d("BP2-BP3") / mean(c(d("BP2-OUT"), d("BP3-OUT"))),
             rnd_TP2_BP3 = d("TP2-BP3") / mean(c(d("TP2-OUT"), d("BP3-OUT"))))
}

bg <- wide[in_peak == FALSE]
block_null <- function(k) {
  b <- copy(bg)
  b[, block := (seq_len(.N) - 1) %/% k, by = chr]
  b[, if (.N == k) summarise_windows(.SD), by = .(chr, block)]
}

res <- rbindlist(lapply(c(2, 3, 6), function(pk) {
  p <- peaks[peak == pk]
  w <- wide[chr == p$chr & start <= p$end & end >= p$start]
  obs <- summarise_windows(w)
  null <- block_null(nrow(w))
  stats <- names(obs)
  data.table(peak = pk, statistic = stats, windows_5kb = nrow(w),
             peak_value = unlist(obs),
             genome_median = sapply(stats, function(s) median(null[[s]][is.finite(null[[s]])])),
             percentile = sapply(stats, function(s) { x <- null[[s]][is.finite(null[[s]])]; 100 * mean(x < obs[[s]]) }),
             null_blocks = nrow(null))
}))

fwrite(res, out_file, sep = "\t")
print(res, digits = 3)
