"""panel_lib.py -- model extension, Stage 1: the summary-statistic panel, shared by observed and simulated data.

Every function takes a genotype matrix + positions + callable layout and nothing data-source specific, so the
same code is applied unchanged to simulations (pass f = 1, callable all True, pol all +1, no missing data).

Inputs per segment (one 1-Mb mask block for observed data; one simulated callable-layout chunk for simulations)
  G         (n_snps, n_ind) int8 ALT dosage 0/1/2, -1 = missing; biallelic SNVs, pos-sorted
  pos       (n_snps,) int64 positions (physical; layout units when f = 1)
  f         callable fraction of the segment (physical distance x f = callable-layout distance; ld23_lib)
  L         callable bp of the segment (denominator of pi and dxy)
  callable_ (n_snps,) bool  SNP lies in the callable set that L counts (diversity, divergence, Fst, delta)
  pol       (n_snps,) int8  +1 derived allele = ALT, -1 derived = REF, 0 unpolarized (drift-test / simulation 24a polarization)
  samples   {name: column-index array}; pairs [(bp_name, tp_name)]

Statistics
  pi        sum over callable SNPs of 2a(n-a)/(n(n-1)) (a ALT copies among n called haplotypes, n >= 2) / L
  dxy       sum over callable SNPs (n1, n2 >= 2) of p1(1-p2) + p2(1-p1) / L
  fst       Hudson Fst, ratio of sums (Bhatia et al. 2013): sum[(p1-p2)^2 - p1(1-p1)/(n1-1) - p2(1-p2)/(n2-1)]
            / sum[p1(1-p2) + p2(1-p1)] over the same SNPs
  delta     f_BP(derived) - f_TP(derived) at callable, polarized SNPs with complete data in both samples (24a rule:
            every individual called); counts at delta >= 0.5, 0.7, 0.78 (24a threshold), 0.8; n_sites = SNPs tested
  ld        ld23_lib.filter_snps + ld23_lib.binned_r2 UNCHANGED (MAF >= 0.10, call >= 0.90 within the sample,
            2-kb bins to 100 kb in layout units) on all SNPs of the segment, as simulation/23a
  ld_far    pooled mean r^2 over 20-100 kb = sum(sum_r2[bins 10..49]) / sum(n_pairs[bins 10..49])
  interchrom r^2 between SNPs on different chromosomes (ld23 filters, >= MIN_JOINT jointly called), SNPs drawn
            uniformly at random, K per chromosome; perm = same with individuals permuted within the second
            chromosome (one permutation per chromosome) -- the 1/n sampling floor; excess = obs - perm
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SIM_DIR = os.environ.get("LD23_DIR", os.path.join(HERE, "..", "..", "..", "simulation"))   # code/11_Drift_vs_Selection/simulation
sys.path.insert(0, SIM_DIR)
import ld23_lib as L  # noqa: E402  (simulation/ld23_lib.py, unchanged)

THRESH = (0.5, 0.7, 0.78, 0.8)
EPS = 1e-9
FAR = slice(10, 50)          # 20-100 kb
K_INTER = 400                # SNPs drawn per chromosome for the interchromosomal statistic
SEED_INTER = 20261108


def _counts(Gs):
    M = Gs >= 0
    X = np.where(M, Gs, 0).astype(np.int64)
    return X.sum(1), 2 * M.sum(1), M


def ld_keep(Gs):
    """ld23 filter for one sample's genotype block; returns (keep, X int8, M int8)."""
    M = (Gs >= 0).astype(np.int8)
    X = np.where(Gs >= 0, Gs, 0).astype(np.int8)
    return L.filter_snps(X, M), X, M


def segment_stats(G, pos, f, L_callable, callable_, pol, samples, pairs):
    """Per-segment sufficient statistics. Returns {(kind, name): value} with arrays for LD bins."""
    out = {}
    cache = {}
    for name, idx in samples.items():
        Gs = G[:, idx]
        a, n, _ = _counts(Gs)
        cache[name] = (a, n, len(idx))
        c = callable_ & (n >= 2)
        with np.errstate(invalid="ignore", divide="ignore"):
            h = np.where(c, 2.0 * a * (n - a) / (n * (n - 1.0)), 0.0)
        out[("S_pi", name)] = float(h.sum())
        out[("n_segsites", name)] = int(np.sum(c & (a > 0) & (a < n)))
        keep, X, M = ld_keep(Gs)
        if keep.sum() >= 2:
            cnt, tot = L.binned_r2(pos[keep], X[keep], M[keep], f)
        else:
            cnt, tot = np.zeros(L.NBIN, np.int64), np.zeros(L.NBIN)
        out[("ld_cnt", name)] = cnt
        out[("ld_tot", name)] = tot
        out[("n_ld_snps", name)] = int(keep.sum())
    for bp, tp in pairs:
        a1, n1, k1 = cache[bp]
        a2, n2, k2 = cache[tp]
        key = f"{bp}|{tp}"
        c = callable_ & (n1 >= 2) & (n2 >= 2)
        with np.errstate(invalid="ignore", divide="ignore"):
            p1, p2 = a1 / n1, a2 / n2
            dxy = p1 * (1 - p2) + p2 * (1 - p1)
            num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1.0) - p2 * (1 - p2) / (n2 - 1.0)
        out[("S_dxy", key)] = float(np.where(c, dxy, 0.0).sum())
        out[("S_fst_num", key)] = float(np.where(c, num, 0.0).sum())
        out[("S_fst_den", key)] = float(np.where(c, dxy, 0.0).sum())
        full = (n1 == 2 * k1) & (n2 == 2 * k2) & (pol != 0)
        for tag, sel in (("", callable_ & full), ("_24arule", full)):
            d1 = np.where(pol == 1, a1, n1 - a1)[sel]
            d2 = np.where(pol == 1, a2, n2 - a2)[sel]
            delta = d1 / (2.0 * k1) - d2 / (2.0 * k2)
            out[("delta_n_sites" + tag, key)] = int(sel.sum())
            for t in THRESH:
                out[(f"delta_ge_{t:g}" + tag, key)] = int(np.sum(delta >= t - EPS))
    out[("L", "")] = float(L_callable)
    return out


def draw_interchrom(rows_by_chrom, seed=SEED_INTER, k=K_INTER):
    """rows_by_chrom: {chrom_number: (X int8 (m, n_ind), M int8 (m, n_ind))} of ld23-filtered SNPs of one sample
    on that chromosome (all mask segments, in position order). Returns {chrom_number: (X, M, Xperm, Mperm)} with
    k SNPs drawn without replacement (rng = default_rng((seed, chrom_number))) and one individual permutation per
    chromosome drawn from the same rng after the SNP draw."""
    out = {}
    for c in sorted(rows_by_chrom):
        X, M = rows_by_chrom[c]
        rng = np.random.default_rng((seed, int(c)))
        m = X.shape[0]
        pick = np.sort(rng.choice(m, min(k, m), replace=False)) if m else np.zeros(0, int)
        perm = rng.permutation(X.shape[1])
        Xs, Ms = X[pick].astype(float), M[pick].astype(float)
        out[c] = (Xs, Ms, Xs[:, perm], Ms[:, perm])
    return out


def _cross(XA, MA, XB, MB):
    # as 23g_unlinked_floor.cross, returning (sum r2, n pairs)
    n = MA @ MB.T; sx = XA @ MB.T; sy = MA @ XB.T
    sxx = (XA ** 2) @ MB.T; syy = MA @ (XB ** 2).T; sxy = XA @ XB.T
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = sxy / n - sx * sy / n ** 2; vx = sxx / n - (sx / n) ** 2; vy = syy / n - (sy / n) ** 2
        r2 = cov ** 2 / (vx * vy)
    ok = (n >= L.MIN_JOINT) & np.isfinite(r2) & (vx > 0) & (vy > 0)
    return float(r2[ok].sum()), int(ok.sum())


def interchrom_matrices(drawn):
    """drawn from draw_interchrom. Returns chrom list and symmetric (C, C) matrices S, N, Sp, Np (zero diagonal):
    sums of r^2 and pair counts between chromosomes, observed and permuted (A unpermuted vs B permuted, A < B)."""
    cs = sorted(drawn)
    C = len(cs)
    S, N, Sp, Np = (np.zeros((C, C)) for _ in range(4))
    for i in range(C):
        XA, MA, _, _ = drawn[cs[i]]
        for j in range(i + 1, C):
            XB, MB, XBp, MBp = drawn[cs[j]]
            s, n = _cross(XA, MA, XB, MB)
            sp, nq = _cross(XA, MA, XBp, MBp)
            S[i, j] = S[j, i] = s; N[i, j] = N[j, i] = n
            Sp[i, j] = Sp[j, i] = sp; Np[i, j] = Np[j, i] = nq
    return cs, S, N, Sp, Np
