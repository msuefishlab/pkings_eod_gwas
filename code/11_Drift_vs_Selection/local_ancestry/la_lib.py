"""la_lib.py -- shared code for the local-ancestry analysis. Not a script.

Two-state haploid HMM (BAM vs TP ancestry) and three-state diploid HMM (0/1/2 BAM copies) on
ancestry-informative markers, with a pulse-admixture transition: between markers d bp apart the
ancestry is redrawn with probability 1 - exp(-lam*d), BAM with probability m; lam = T * r, fitted at
r = 1.4 cM/Mb (R_FIT per bp). Emissions: empirical-Bayes panel frequencies (beta prior from 03_markers.py)
with error EPS.

Grid likelihoods are computed per chromosome, vectorized over (T, m, haplotype or fish), with scaled forward
recursions. Posteriors use per-haplotype (or per-fish) parameters.
"""
import numpy as np

EPS = 0.01
R_FIT = 1.4e-8                      # Morgans per bp (1.4 cM/Mb)
R_BRACKET = (1.0, 1.4, 3.0)         # cM/Mb
T_GRID = np.geomspace(1, 20000, 41)
M_GRID = np.array([0.001, 0.003] + [i / 100 for i in range(1, 100)] + [0.997, 0.999])
INV = ("chr4", 40_880_000, 42_770_000)


def load_markers(out, panel):
    z = np.load(f"{out}/markers_{panel}.npz", allow_pickle=False)
    d = {k: z[k] for k in z.files}
    o = d["orient"].astype(bool)
    (aB, bB), (aT, bT) = d["prior_B"], d["prior_T"]
    # empirical-Bayes frequency of the Bambomo-high allele, then back to allele 1
    qB = (np.where(o, d["cB"], d["nB"] - d["cB"]) + aB) / (d["nB"] + aB + bB)
    qT = (np.where(o, d["cT"], d["nT"] - d["cT"]) + aT) / (d["nT"] + aT + bT)
    pB = np.where(o, qB, 1 - qB)
    pT = np.where(o, qT, 1 - qT)
    d["pBe"] = pB * (1 - EPS) + (1 - pB) * EPS
    d["pTe"] = pT * (1 - EPS) + (1 - pT) * EPS
    return d


def excluded_regions():
    excl = [INV]
    for line in open("input_data/10_Introgression/peaks.bed"):
        f = line.split()
        excl.append((f"chr{f[0]}", int(f[1]) - 250_000, int(f[2]) + 250_000))
    return excl


def segments(chrom_num, pos, excl):
    """Segment id per marker: a new segment starts on a chromosome change or when an excluded interval lies
    between two consecutive markers."""
    seg = np.zeros(pos.size, int)
    s = 0
    for i in range(1, pos.size):
        if chrom_num[i] != chrom_num[i - 1]:
            s += 1
        else:
            ch = f"chr{chrom_num[i]}"
            if any(x == ch and pos[i - 1] < a and pos[i] > b for x, a, b in excl):
                s += 1
        seg[i] = s
    return seg


# ---------------------------------------------------------------- haploid, grid likelihood
def hap_loglik_grid(obs, pos, pBe, pTe, t_grid=T_GRID, m_grid=M_GRID):
    """obs: (n, H) 0/1 alleles on one chromosome. Returns log-likelihood (len(t_grid), len(m_grid), H)."""
    n, H = obs.shape
    lam = (t_grid * R_FIT)[:, None, None]
    m = m_grid[None, :, None]
    eB = np.where(obs == 1, pBe[:, None], 1 - pBe[:, None])
    eT = np.where(obs == 1, pTe[:, None], 1 - pTe[:, None])
    shape = (t_grid.size, m_grid.size, H)
    aB = np.broadcast_to(m * eB[0], shape).copy()
    aT = np.broadcast_to((1 - m) * eT[0], shape).copy()
    s = aB + aT
    ll = np.log(s); aB /= s; aT /= s
    c = np.ones(shape)
    d = np.diff(pos).astype(float)
    for j in range(1, n):
        e = np.exp(-lam * d[j - 1])
        rm = (1 - e) * m
        aB *= e; aB += rm; aB *= eB[j]
        aT *= e; aT += (1 - e) - rm; aT *= eT[j]
        s = aB + aT
        aB /= s; aT /= s; c *= s
        if j % 25 == 0:
            ll += np.log(c); c[:] = 1
    ll += np.log(c)
    return ll


# ---------------------------------------------------------------- diploid, grid likelihood
def dip_emissions(geno, pBe, pTe):
    """geno: (n, F) 0/1/2, -1 missing. Returns E (n, F, 3) for k = 0, 1, 2 BAM copies."""
    pb, pt = pBe[:, None], pTe[:, None]
    E = np.empty(geno.shape + (3,))
    for g, (f0, f1, f2) in enumerate((((1 - pt) ** 2, (1 - pb) * (1 - pt), (1 - pb) ** 2),
                                      (2 * pt * (1 - pt), pb * (1 - pt) + (1 - pb) * pt, 2 * pb * (1 - pb)),
                                      (pt ** 2, pb * pt, pb ** 2))):
        msk = geno == g
        E[..., 0] = np.where(msk, f0, E[..., 0]) if g else np.where(msk, f0, 1.0)
        E[..., 1] = np.where(msk, f1, E[..., 1]) if g else np.where(msk, f1, 1.0)
        E[..., 2] = np.where(msk, f2, E[..., 2]) if g else np.where(msk, f2, 1.0)
    return E


def _dip_trans(e, m):
    qBB = e + (1 - e) * m
    qTB = (1 - e) * m
    P = {}
    P[0, 0] = (1 - qTB) ** 2; P[0, 1] = 2 * qTB * (1 - qTB); P[0, 2] = qTB ** 2
    P[2, 2] = qBB ** 2; P[2, 1] = 2 * qBB * (1 - qBB); P[2, 0] = (1 - qBB) ** 2
    P[1, 2] = qBB * qTB; P[1, 0] = (1 - qBB) * (1 - qTB); P[1, 1] = qBB * (1 - qTB) + (1 - qBB) * qTB
    return P


def dip_loglik_grid(geno, pos, pBe, pTe, t_grid=T_GRID, m_grid=M_GRID):
    """Returns log-likelihood (len(t_grid), len(m_grid), F)."""
    n, F = geno.shape
    E = dip_emissions(geno, pBe, pTe)
    lam = (t_grid * R_FIT)[:, None, None]
    m = m_grid[None, :, None]
    shape = (t_grid.size, m_grid.size, F)
    a = [np.broadcast_to(w * E[0, :, k], shape).copy() for k, w in enumerate(((1 - m) ** 2, 2 * m * (1 - m), m ** 2))]
    s = a[0] + a[1] + a[2]
    ll = np.log(s); a = [x / s for x in a]
    c = np.ones(shape)
    d = np.diff(pos).astype(float)
    for j in range(1, n):
        P = _dip_trans(np.exp(-lam * d[j - 1]), m)
        nxt = [(a[0] * P[0, k] + a[1] * P[1, k] + a[2] * P[2, k]) * E[j, :, k] for k in range(3)]
        s = nxt[0] + nxt[1] + nxt[2]
        a = [x / s for x in nxt]; c *= s
        if j % 25 == 0:
            ll += np.log(c); c[:] = 1
    ll += np.log(c)
    return ll


# ---------------------------------------------------------------- posteriors with per-column parameters
def hap_posterior(obs, pos, pBe, pTe, lam, m):
    """obs (n, H); lam, m (H,). Returns P(BAM) (n, H)."""
    n, H = obs.shape
    eB = np.where(obs == 1, pBe[:, None], 1 - pBe[:, None])
    eT = np.where(obs == 1, pTe[:, None], 1 - pTe[:, None])
    d = np.diff(pos).astype(float)
    fB = np.empty((n, H)); fT = np.empty((n, H))
    aB, aT = m * eB[0], (1 - m) * eT[0]
    s = aB + aT; fB[0], fT[0] = aB / s, aT / s
    E = np.exp(-lam[None, :] * d[:, None]) if n > 1 else np.zeros((0, H))
    for j in range(1, n):
        e = E[j - 1]
        aB = (e * fB[j - 1] + (1 - e) * m) * eB[j]
        aT = (e * fT[j - 1] + (1 - e) * (1 - m)) * eT[j]
        s = aB + aT; fB[j], fT[j] = aB / s, aT / s
    bB = np.ones(H); bT = np.ones(H)
    post = np.empty((n, H))
    post[n - 1] = fB[n - 1]
    for j in range(n - 2, -1, -1):
        e = E[j]
        xB, xT = bB * eB[j + 1], bT * eT[j + 1]
        mix = (1 - e) * (m * xB + (1 - m) * xT)
        nB, nT = e * xB + mix, e * xT + mix
        s = nB + nT; bB, bT = nB / s, nT / s
        g = fB[j] * bB; h = fT[j] * bT
        post[j] = g / (g + h)
    return post


def dip_posterior(geno, pos, pBe, pTe, lam, m):
    """geno (n, F); lam, m (F,). Returns state posteriors (n, F, 3)."""
    n, F = geno.shape
    E = dip_emissions(geno, pBe, pTe)
    d = np.diff(pos).astype(float)
    f = np.empty((n, F, 3))
    a = np.stack([(1 - m) ** 2, 2 * m * (1 - m), m ** 2], -1) * E[0]
    f[0] = a / a.sum(-1, keepdims=True)
    Ps = []
    for j in range(1, n):
        P = _dip_trans(np.exp(-lam * d[j - 1]), m)
        M = np.empty((F, 3, 3))
        for (x, y), v in P.items():
            M[:, x, y] = v
        Ps.append(M)
        a = np.einsum("fi,fij->fj", f[j - 1], M) * E[j]
        f[j] = a / a.sum(-1, keepdims=True)
    b = np.ones((F, 3))
    post = np.empty((n, F, 3))
    post[n - 1] = f[n - 1]
    for j in range(n - 2, -1, -1):
        b = np.einsum("fij,fj->fi", Ps[j], b * E[j + 1])
        b /= b.sum(-1, keepdims=True)
        g = f[j] * b
        post[j] = g / g.sum(-1, keepdims=True)
    return post


# ---------------------------------------------------------------- fitting helpers
def fit_from_surfaces(LLc, B=1000, seed=20261103):
    """LLc: (n_chrom, G1, G2, F) per-chromosome log-likelihood surfaces. Returns per fish dict with ML indices
    and bootstrap (chromosome resampling) 2.5/97.5 percentiles of T and m."""
    nc, G1, G2, F = LLc.shape
    tot = LLc.sum(0)
    rng = np.random.default_rng(seed)
    W = np.stack([np.bincount(rng.integers(0, nc, nc), minlength=nc) for _ in range(B)]).astype(float)
    res = []
    for f in range(F):
        i, k = np.unravel_index(np.argmax(tot[..., f]), (G1, G2))
        boot = W @ LLc[..., f].reshape(nc, -1)
        bi, bk = np.unravel_index(np.argmax(boot, 1), (G1, G2))
        res.append(dict(iT=i, im=k, T=T_GRID[i], m=M_GRID[k], ll=tot[i, k, f],
                        T_lo=np.percentile(T_GRID[bi], 2.5), T_hi=np.percentile(T_GRID[bi], 97.5),
                        m_lo=np.percentile(M_GRID[bk], 2.5), m_hi=np.percentile(M_GRID[bk], 97.5),
                        T_edge=i in (0, G1 - 1), m_edge=k in (0, G2 - 1)))
    return res


def call_tracts(post, pos, seg, thr=0.5):
    """post (n,) P(BAM) for one haplotype; pos, seg (n,). Returns list of (seg, start, end) tracts with edges at
    the midpoints to the flanking markers of the same segment (segment ends at the end markers)."""
    out = []
    n = pos.size
    on = post >= thr
    i = 0
    while i < n:
        if not on[i]:
            i += 1; continue
        j = i
        while j + 1 < n and on[j + 1] and seg[j + 1] == seg[i]:
            j += 1
        left = pos[i] if (i == 0 or seg[i - 1] != seg[i]) else (pos[i - 1] + pos[i]) / 2
        right = pos[j] if (j == n - 1 or seg[j + 1] != seg[j]) else (pos[j] + pos[j + 1]) / 2
        out.append((seg[i], left, right))
        i = j + 1
    return out


def analysed_length(pos, seg):
    return float(sum(pos[seg == s].max() - pos[seg == s].min() for s in np.unique(seg)))


# ---------------------------------------------------------------- orchestration shared by 04-06
def _surface_worker(args):
    mode, obs, pos, pBe, pTe = args
    f = hap_loglik_grid if mode == "hap" else dip_loglik_grid
    return f(obs, pos, pBe, pTe)


def grid_surfaces(mode, per_chrom, workers=12):
    """per_chrom: list of (obs, pos, pBe, pTe), one per chromosome. Returns (n_chrom, G1, G2, columns)."""
    import multiprocessing as mp
    with mp.get_context("fork").Pool(workers) as pool:
        res = pool.map(_surface_worker, [(mode,) + tuple(x) for x in per_chrom], chunksize=1)
    return np.stack(res)


def per_chrom_split(chrom, pos, cols, pBe, pTe, mask):
    out = []
    for c in range(1, 26):
        s = (chrom == c) & mask
        out.append((cols[s], pos[s], pBe[s], pTe[s]))
    return out


def summarise_fish(post_h, pos, chrom, seg, T, m, r_fit_cm=1.4):
    """post_h: (n, 2) P(BAM) on the fish's two haplotypes (genome-wide markers). Returns summary dict and the
    tract list [(hap, chrom, start, end)]."""
    tracts = []
    for h in range(post_h.shape[1]):
        for s, a, b in call_tracts(post_h[:, h], pos, seg):
            tracts.append((h, int(chrom[seg == s][0]), a, b))
    L = np.array([b - a for _, _, a, b in tracts]) if tracts else np.zeros(0)
    alen = analysed_length(pos, seg) * post_h.shape[1]
    frac_tract = L.sum() / alen
    Lbar = L.mean() if L.size else np.nan
    r = r_fit_cm * 1e-8
    T_tract = 1 / ((1 - m) * r * Lbar) if L.size and m < 1 else np.nan
    F_long = float(L[L >= 2_000_000].sum() / alen) if L.size else 0.0
    return dict(frac_post=float(post_h.mean()), n_tracts=int(L.size), F_long=F_long,
                med_tract=float(np.median(L)) if L.size else np.nan, mean_tract=float(Lbar) if L.size else np.nan,
                tot_tract=float(L.sum()), max_tract=float(L.max()) if L.size else 0.0,
                frac_tract=float(frac_tract), T_from_tracts=float(T_tract)), tracts


# ---------------------------------------------------------------- statistics and simulation matching
LONG_BP = 2_000_000
SD_FLOOR = {"q_m": 0.01, "m_hat": 0.01, "logT": 0.05, "F_long": 0.01, "m_hat_dip": 0.01, "logT_dip": 0.05}
MATCH_KEYS = {"phased": ("q_m", "m_hat", "logT", "F_long"), "diploid": ("q_m", "m_hat_dip", "logT_dip")}


def projection_q(x, pBe, pTe):
    """x: (n, k) allele dosages per haplotype (0/1) or genotype/2 (0..1); -1 or nan missing. T-free BAM fraction."""
    x = np.asarray(x, float)
    ok = (x >= 0) & np.isfinite(x)
    dif = (pBe - pTe)[:, None] * np.ones_like(x)
    num = np.where(ok, (x - pTe[:, None]) * dif, 0).sum()
    den = np.where(ok, dif * dif, 0).sum()
    return float(num / den)


def long_fraction(tracts, alen_total, min_bp=LONG_BP):
    L = np.array([b - a for _, _, a, b in tracts]) if tracts else np.zeros(0)
    return float(L[L >= min_bp].sum() / alen_total)


def cell_stats(sim_rows, keys):
    """sim_rows: list of dicts with 'T', 'm' and the statistics. Returns {(T, m): {key: values array}}."""
    cells = {}
    for r in sim_rows:
        cells.setdefault((r["T"], r["m"]), []).append(r)
    return {c: {k: np.array([float(r[k]) for r in v]) for k in keys} for c, v in cells.items()}


def match(x, cells, keys, exclude=None):
    """x: dict of the fish's statistics. cells from cell_stats. exclude: (cell, index) to leave out (LOO).
    Returns (best_cell, compatible_cells, chi2 by cell)."""
    chi = {}
    compat = []
    for c, st in cells.items():
        zz = []
        for k in keys:
            v = st[k]
            if exclude is not None and exclude[0] == c:
                v = np.delete(v, exclude[1])
            sd = max(v.std(ddof=1) if v.size > 1 else 0.0, SD_FLOOR[k])
            zz.append((float(x[k]) - v.mean()) / sd)
        zz = np.array(zz)
        chi[c] = float((zz ** 2).sum())
        if np.all(np.abs(zz) <= 3):
            compat.append(c)
    best = min(chi, key=chi.get)
    return best, compat, chi
