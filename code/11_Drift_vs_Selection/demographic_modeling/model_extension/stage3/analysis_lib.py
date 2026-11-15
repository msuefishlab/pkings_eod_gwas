"""analysis_lib.py -- Stage 3 fitting and scoring (run locally on the compact step outputs).

r/mu fit: the 23e procedure, with 23e's own interp_family / fit functions imported unchanged: observed curve = pooled
sum r^2 / sum n_pairs over the chromosomes of the fit set; simulated family = pooled curve per grid value; linear
interpolation in log r/mu on 2,000 points; SSE minimum; 1,000 chromosome-bootstrap refits (resamples = default_rng(
20261030), one integers(0, n, n) per resample, fresh rng per (population, set, observed row) -- the Stage 1
ext_aggregate draw, so the band and the refits use the same resamples); 95% percentile interval; edge = argmin at the
first or last fine-grid point; 'shape' = every curve divided by its first bin.
"""
import csv, importlib.util, os
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("e23", os.path.join(HERE, "..", "..", "..", "simulation",
                                                                    "23e_fit_r_over_mu.py"))
E23 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(E23)
NB, NBOOT, SEED = 50, 1000, 20261030
GRID = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0]
FAR = slice(10, 50)


def read_split(fn):
    return {r["chrom"]: r["set"] for r in csv.DictReader(open(fn), delimiter="\t")}


def chrom_sums(npz, pop, chroms):
    z = np.load(npz)
    N = np.vstack([z[f"{c}|{pop}|ld_cnt"] for c in chroms]); T = np.vstack([z[f"{c}|{pop}|ld_tot"] for c in chroms])
    return N, T


def set_chroms(split, s):
    cs = sorted(split) if s == "genome" else sorted(c for c in split if split[c] == s)
    return cs


def boot_idx(n):
    rng = np.random.default_rng(SEED)
    return np.vstack([rng.integers(0, n, n) for _ in range(NBOOT)])


def obs_curve_boot(N, T):
    obs = T.sum(0) / N.sum(0)
    idx = boot_idx(N.shape[0])
    bt = np.vstack([T[i].sum(0) / N[i].sum(0) for i in idx])
    lo, hi = np.percentile(bt, [2.5, 97.5], axis=0)
    return obs, bt, lo, hi


def family(fs, m, fam, t, pop):
    """{grid value: pooled curve} for model m, family fam ('cal' | 'ins' | 'ho' (= sets 1-5 pooled) | 'ho_sK')"""
    out = {}
    for k in range(1, 13):
        if fam == "ho":
            c = sum(fs[f"{m}|ho_s{s}|{k}|{t}|{pop}|cnt"] for s in range(1, 6))
            tt = sum(fs[f"{m}|ho_s{s}|{k}|{t}|{pop}|tot"] for s in range(1, 6))
        else:
            c, tt = fs[f"{m}|{fam}|{k}|{t}|{pop}|cnt"], fs[f"{m}|{fam}|{k}|{t}|{pop}|tot"]
        out[GRID[k - 1]] = tt / c
    return out


def family_far(fs, m, t, pop, bins=FAR):
    out = {}
    for k in range(1, 13):
        c = sum(fs[f"{m}|ho_s{s}|{k}|{t}|{pop}|cnt"] for s in range(1, 6))
        tt = sum(fs[f"{m}|ho_s{s}|{k}|{t}|{pop}|tot"] for s in range(1, 6))
        out[GRID[k - 1]] = tt[bins].sum() / c[bins].sum()
    return out


def interp_at(grid, r):
    """23e interpolation (linear in log r/mu between grid points) of a {r: curve-or-scalar} family at r (clipped)"""
    ks = np.array(sorted(grid)); C = np.array([grid[k] for k in ks], float)
    lr = np.log(np.clip(r, ks[0], ks[-1])); lk = np.log(ks)
    j = int(np.clip(np.searchsorted(lk, lr, "right") - 1, 0, len(ks) - 2))
    w = (lr - lk[j]) / (lk[j + 1] - lk[j])
    return (1 - w) * C[j] + w * C[j + 1]


def fit23e(N, T, grid, shape=False, boots=True):
    """returns dict(r_hat, lo95, hi95, edge, frac_in_band(own set band), curve_at_fit, obs, lo, hi)"""
    obs, bt, lo, hi = obs_curve_boot(N, T)
    fine, F = E23.interp_family(grid)
    if shape:
        F = F / F[:, :1]; o = obs / obs[0]; bk = bt / bt[:, :1]
    else:
        o, bk = obs, bt
    k_hat, i_hat = E23.fit(o, fine, F)
    if boots:
        kb = np.array([E23.fit(bk[j], fine, F)[0] for j in range(NBOOT)])
        klo, khi = np.percentile(kb, [2.5, 97.5])
    else:
        klo = khi = np.nan
    blo, bhi = np.percentile(bk, [2.5, 97.5], axis=0)
    inside = float(np.mean((F[i_hat] >= blo) & (F[i_hat] <= bhi)))
    return dict(r_hat=float(k_hat), lo95=float(klo), hi95=float(khi), edge=int(i_hat in (0, len(fine) - 1)),
                frac_in_band_own=inside, curve=F[i_hat], obs=o, lo=blo, hi=bhi,
                sse=float(((F[i_hat] - o) ** 2).sum()))


def bins_in_band(curve, lo, hi):
    inb = (curve >= lo) & (curve <= hi)
    return float(np.mean(inb)), int(np.sum(curve > hi)), int(np.sum(curve < lo))


def clopper_pearson(k, n, a=0.05):
    from scipy.stats import beta
    lo = 0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k)
    return float(lo), float(hi)


def read_panel(fn):
    out = {}
    for r in csv.DictReader(open(fn), delimiter="\t"):
        out[(r["statistic"], r["sample"], r["set"])] = (float(r["value"]), float(r["lo"]), float(r["hi"]))
    return out


def read_curves(fn):
    out = defaultdict(dict)
    for r in csv.DictReader(open(fn), delimiter="\t"):
        if r["kind"] == "observed":
            out[(r["pop"], r["set"])][int(r["bin"])] = (float(r["r2"]), float(r["lo"]), float(r["hi"]),
                                                       float(r["n_pairs"]))
    res = {}
    for k, d in out.items():
        a = np.array([d[b] for b in range(NB)])
        res[k] = dict(r2=a[:, 0], lo=a[:, 1], hi=a[:, 2], n=a[:, 3])
    return res
