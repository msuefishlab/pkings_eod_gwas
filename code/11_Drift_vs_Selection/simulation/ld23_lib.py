"""ld23_lib.py -- the one LD-decay statistic shared by 23a (observed) and 23b (simulated).

Genotype r^2: the squared Pearson correlation of 0/1/2 dosages over the individuals called at both SNPs,
for SNP pairs closer than MAXD in callable-layout distance, summed per BIN-wide distance bin. Callable-layout
distance is physical distance x the callable fraction of the 1-Mb block (the simulation lays the callable genome
out contiguously, so its distances are already in these units: f = 1).

Both callers keep SNPs with within-sample MAF >= MIN_MAF and call rate >= MIN_CALL. Simulated genotypes have no
missing data, so for them the jointly-called formula below reduces to the ordinary correlation.
"""
import numpy as np

MAXD = 100_000          # bp, callable-layout units
BIN = 2_000
NBIN = MAXD // BIN      # 50 bins
MIN_MAF = 0.10
MIN_CALL = 0.90
MIN_JOINT = 10          # a pair needs this many jointly called individuals
CHUNK = 400


def filter_snps(X, M):
    """X: (n_snps, n_ind) dosages with 0 where missing; M: (n_snps, n_ind) 1 if called.
    Returns the boolean keep mask (call rate and within-sample MAF)."""
    ncall = M.sum(1)
    call_rate = ncall / M.shape[1]
    with np.errstate(invalid="ignore", divide="ignore"):
        p = X.sum(1) / (2.0 * ncall)
    maf = np.minimum(p, 1 - p)
    return (call_rate >= MIN_CALL) & (maf >= MIN_MAF)


def binned_r2(pos, X, M, f=1.0):
    """pos: sorted positions (n,); X, M as in filter_snps, already filtered; f: callable fraction.
    Returns (n_pairs[NBIN], sum_r2[NBIN])."""
    n = len(pos)
    cnt = np.zeros(NBIN, np.int64)
    tot = np.zeros(NBIN, np.float64)
    if n < 2:
        return cnt, tot
    X = X.astype(np.float64)
    M = M.astype(np.float64)
    X2 = X * X
    maxphys = MAXD / f
    for i0 in range(0, n, CHUNK):
        i1 = min(n, i0 + CHUNK)
        j1 = int(np.searchsorted(pos, pos[i1 - 1] + maxphys, "right"))
        A, B = slice(i0, i1), slice(i0, j1)
        nij = M[A] @ M[B].T
        sx, sy = X[A] @ M[B].T, M[A] @ X[B].T
        sxx, syy = X2[A] @ M[B].T, M[A] @ X2[B].T
        sxy = X[A] @ X[B].T
        with np.errstate(invalid="ignore", divide="ignore"):
            cov = sxy / nij - sx * sy / nij ** 2
            vx = sxx / nij - (sx / nij) ** 2
            vy = syy / nij - (sy / nij) ** 2
            r2 = cov * cov / (vx * vy)
        gi = np.arange(i0, i1)[:, None]
        gj = np.arange(i0, j1)[None, :]
        d = (pos[B][None, :] - pos[A][:, None]) * f
        ok = (gj > gi) & (d < MAXD) & (nij >= MIN_JOINT) & np.isfinite(r2) & (vx > 0) & (vy > 0)
        b = (d[ok] // BIN).astype(int)
        cnt += np.bincount(b, minlength=NBIN)[:NBIN]
        tot += np.bincount(b, weights=r2[ok], minlength=NBIN)[:NBIN]
    return cnt, tot
