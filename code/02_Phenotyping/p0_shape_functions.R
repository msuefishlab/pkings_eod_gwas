# P0 shape detectors, shared by:
#   - p0_phenotype_assessment.Rmd    (the by-eye-calibrated two-axis classifier)
#   - p0_statistical_detection.Rmd   (the label-free per-fish detection layer)
#
# Both must score P0 windows with the IDENTICAL detector, so the primitives live
# here. Two families:
#
#   1. make_bank / score_matrix / snr_matrix  -- verbatim originals, bank sized to
#      the input window width. Used as-is by p0_phenotype_assessment.Rmd.
#   2. p0_score_one / slide_scores / p0_scores_canonical  -- a canonical-length
#      wrapper (default L = 128) so an observed P0 window and every sliding
#      baseline sub-window get scored the same way regardless of native length.
#      This is what the empirical-null layer needs: the observed window and the
#      null strip must be treated identically or the p-value is meaningless.
#
# The presence signals are amplitude-invariant by construction: `shape` is a
# max Pearson correlation (scale-free) and `snr` divides the dip depth by the
# window's own MAD. Combined with a per-fish null, the population size confound
# in absolute P0 amplitude cannot leak into the presence call.

# ---- 1. Original detectors (bank sized to input width) ---------------------

# Bank of small-dip templates spanning real P0 positions (they sit right up
# against the P1 onset, ~0.80-0.98 of the window) and widths. A single fixed
# template misses edge-positioned / broad dips; we take the best fit over the bank.
# Positions run to 0.98: several real wobbles (APA_207/309/6681, BAM_6550) carry a
# sharp P0 dip jammed right against the P1 onset at 0.96-0.98 of the window, which a
# bank ending at 0.94 would miss (deflating their shape score by up to 0.31). The
# near-edge templates only fire on a genuine edge-dip -- they do NOT raise any BP
# fish's shape (BP has no coherent dip to match), so BP specificity is unaffected.
make_bank <- function(n) {
  x <- seq(0, 1, length.out = n)
  bank <- list()
  for (pos in seq(0.80, 0.98, 0.02))
    for (wd in c(0.06, 0.10, 0.15)) {
      t <- -exp(-((x - pos) / wd)^2)
      bank[[length(bank) + 1]] <- t - mean(t)
    }
  bank
}

# Best-matching template per z-scored window (max correlation over the bank).
score_matrix <- function(mat) {
  bank <- make_bank(ncol(mat))
  apply(mat, 1, function(w) {
    w <- w - mean(w)
    if (sd(w) == 0) return(NA_real_)
    max(vapply(bank, function(t) suppressWarnings(cor(w, t)), numeric(1)))
  })
}

# Dip depth relative to the fish's own baseline noise (robust MAD units).
snr_matrix <- function(mat) {
  apply(mat, 1, function(w) {
    n <- length(w)
    base <- w[1:round(0.55 * n)]      # flat left portion of the window
    zone <- w[round(0.6 * n):n]       # P0 zone, right up against the P1 onset
    (median(base) - min(zone)) / (mad(base) + 1e-9)
  })
}

# ---- 2. Canonical-length wrappers for the empirical-null layer -------------

# Resample a numeric vector to length L by linear interpolation.
.resample_to <- function(w, L) {
  if (length(w) == L) return(w)
  approx(seq(0, 1, length.out = length(w)), w, seq(0, 1, length.out = L))$y
}

# Score one window at canonical length L. Resample -> z-score (matching the
# MATLAB p0_window treatment) -> {shape = max corr over bank, snr = dipSNR}.
# Pass a precomputed bank = make_bank(L) to avoid rebuilding it in a loop.
p0_score_one <- function(w, bank, L = 128L) {
  w <- as.numeric(w); w <- w[is.finite(w)]
  if (length(w) < 3) return(c(shape = NA_real_, snr = NA_real_))
  w <- .resample_to(w, L)
  w <- w - mean(w); s <- sd(w); if (s > 0) w <- w / s
  shape <- if (sd(w) == 0) NA_real_
           else max(vapply(bank, function(t) suppressWarnings(cor(w, t)), numeric(1)))
  n <- length(w)
  base <- w[1:round(0.55 * n)]
  zone <- w[round(0.6 * n):n]
  snr <- (median(base) - min(zone)) / (mad(base) + 1e-9)
  c(shape = shape, snr = snr)
}

# Score a matrix of observed windows (rows) at canonical length L.
p0_scores_canonical <- function(mat, L = 128L) {
  bank <- make_bank(L)
  sc <- t(apply(mat, 1, p0_score_one, bank = bank, L = L))
  data.frame(shape = sc[, "shape"], snr = sc[, "snr"])
}

# Slide a length-L window across a baseline strip and score every position.
# Returns the per-fish NULL sample: data.frame(shape, snr) with one row per
# sliding position. `strip` is one row of p0_baseline(_per_sweep).csv (b-cols).
slide_scores <- function(strip, bank, L = 128L, step = 4L) {
  strip <- as.numeric(strip); strip <- strip[is.finite(strip)]
  if (length(strip) < L + 1) return(data.frame(shape = numeric(0), snr = numeric(0)))
  starts <- seq(1L, length(strip) - L + 1L, by = step)
  sc <- t(vapply(starts, function(p) p0_score_one(strip[p:(p + L - 1L)], bank, L),
                 numeric(2)))
  data.frame(shape = sc[, 1], snr = sc[, 2])
}

# iZC1-anchored per-sweep REPRODUCIBILITY of the P0 dip. A = sweeps x points, the
# iZC1-anchored per-sweep P0 windows (rows of p0_windows_per_sweep.csv). A real P0
# recurs across sweeps, so random half-averages of the DIP ZONE (right of the window,
# up against the P1 onset) correlate highly; a BP look-alike's average "dip" is
# non-reproducing noise and the halves disagree (corr ~ 0). Returns the mean
# split-half correlation over `nrep` random splits; needs >= 4 sweeps (else NA).
#
# This is the presence signal that shape/dipSNR (computed on the sweep-AVERAGE) throw
# away: a noise dip can fake a good average shape, but it cannot reproduce across
# sweeps. It is amplitude-invariant (correlation) and, unlike the P1-aligned split-half
# (which is at chance because the P0 is invisible P1-aligned), it operates in the
# iZC1-anchored frame where the P0 power lives. It is also distinct from onset
# phase-locking (tp0_sd): it scores dip SHAPE agreement -- anatomy-confirmed BP scores
# ~0 here while TP scores ~1, and it rescues the anatomy-confirmed wobble APA_6681.
p0_repro_izc1 <- function(A, zone = c(0.60, 1.0), nrep = 200L) {
  A <- A[stats::complete.cases(A), , drop = FALSE]
  ns <- nrow(A); if (ns < 4) return(NA_real_)
  m <- ncol(A); idx <- max(1L, round(zone[1] * m)):round(zone[2] * m)
  Z <- A[, idx, drop = FALSE]
  r <- numeric(nrep)
  for (i in seq_len(nrep)) {
    p <- sample(ns); h1 <- p[1:floor(ns / 2)]; h2 <- p[(floor(ns / 2) + 1):ns]
    r[i] <- suppressWarnings(cor(colMeans(Z[h1, , drop = FALSE]),
                                 colMeans(Z[h2, , drop = FALSE])))
  }
  mean(r, na.rm = TRUE)
}

# ---- 3. Evoked-response detectors on P1-aligned raw traces ------------------
#
# These operate on p0_aligned_per_sweep.csv: A = sweeps x time, ptp-normalized,
# NOT z-scored, on a fixed P1-relative ms axis. Each replicate sweep is a trial of
# the same discharge, so a real P0 is phase-locked (recurs at the same latency) and
# noise is not -- the evoked-response detection setting (ABR/ERP: F_sp, +/- ref).

# Time (ms rel. P1) axis of the aligned export; MUST match the MATLAB constants
# (ALIGN_MS0/ALIGN_MS1/N_align in export_per_sweep_measurements.m).
p0_aligned_axis <- function(n = 448L, ms0 = -3.0, ms1 = 0.5) seq(ms0, ms1, length.out = n)

# Mexican-hat / Ricker DIP template (negative 2nd-derivative-of-Gaussian), the
# wavelet that IS the matched filter for a Gaussian-shaped dip. Centered at mu (ms)
# with scale sig (ms); mean-centered and unit-norm so <trace, w> is a matched-filter
# output invariant to DC and template scaling.
ricker <- function(x, mu, sig) {
  z <- (x - mu) / sig
  t <- -(1 - z^2) * exp(-z^2 / 2)      # negative Ricker = a downward dip
  t <- t - mean(t)
  nn <- sqrt(sum(t^2)); if (nn > 0) t / nn else t
}
ricker_bank <- function(x, mus, sigs) {
  bank <- list()
  for (mu in mus) for (sg in sigs) bank[[length(bank) + 1]] <- ricker(x, mu, sg)
  bank
}
# Max matched-filter projection of a (mean-centered) trace onto a ricker bank.
proj_max <- function(v, bank) {
  vc <- v - mean(v)
  max(vapply(bank, function(w) sum(vc * w), numeric(1)))
}

# Sign-flip (+/- reference) null engine. stat_fun(avg_vector) -> scalar T. Randomly
# flips each sweep's sign before averaging: under H0 (no phase-locked signal) this
# leaves the distribution unchanged; under H1 the signal cancels, so the null T is
# smaller. One-sided empirical p. Robust to autocorrelation (no F-distribution df).
signflip_p <- function(A, stat_fun, nperm = 500L) {
  ns <- nrow(A)
  obs <- stat_fun(colMeans(A))
  if (!is.finite(obs)) return(c(stat = NA_real_, p = NA_real_))
  null <- numeric(nperm)
  for (b in seq_len(nperm)) null[b] <- stat_fun(colMeans(sample(c(-1, 1), ns, TRUE) * A))
  c(stat = obs, p = (1 + sum(null >= obs, na.rm = TRUE)) / (1 + nperm))
}

# F_sp-style energy detector: signal energy of the averaged trace in the P0 window,
# baseline-referenced. Shape-agnostic (any phase-locked deflection counts).
fsp_stat <- function(p0_idx, base_idx) function(avg) {
  a <- avg - mean(avg[base_idx]); sum(a[p0_idx]^2)
}
# Matched-wavelet detector: max ricker-dip projection in the P0 window (shape- and
# landmark-specific). Same +/- null as F_sp.
pm_stat <- function(p0_idx, base_idx, bank) function(avg) {
  a <- avg - mean(avg[base_idx]); proj_max(a[p0_idx], bank)
}

# Split-half reproducibility on the aligned P0 window (score = mean split-half
# correlation; p from a latency-shuffle null). Needs >= 4 sweeps.
splithalf_repro <- function(A, p0_idx, nrep = 200L) {
  ns <- nrow(A); if (ns < 4) return(c(stat = NA_real_, p = NA_real_))
  Z <- A[, p0_idx, drop = FALSE]; m <- length(p0_idx)
  obsr <- nulr <- numeric(nrep)
  for (r in seq_len(nrep)) {
    idx <- sample(ns); h1 <- idx[1:floor(ns / 2)]; h2 <- idx[(floor(ns / 2) + 1):ns]
    a <- colMeans(Z[h1, , drop = FALSE]); b <- colMeans(Z[h2, , drop = FALSE])
    obsr[r] <- suppressWarnings(cor(a, b))
    sh <- sample(2:(m - 1), 1); b2 <- b[c((sh + 1):m, 1:sh)]
    nulr[r] <- suppressWarnings(cor(a, b2))
  }
  rm <- mean(obsr, na.rm = TRUE)
  c(stat = rm, p = (1 + sum(nulr >= rm, na.rm = TRUE)) / (1 + sum(is.finite(nulr))))
}
