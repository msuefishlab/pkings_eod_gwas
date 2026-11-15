#!/usr/bin/env python3
"""23e_fit_r_over_mu.py -- fit r/mu to the observed LD decay under each fitted history (LD-decay test).

Observed curve: mean r^2 per bin pooled over the 200 segments of 23a (sum of sum_r2 / sum of n_pairs).
Simulated curves: the same pooling over the 100 replicates of 23b, one curve per (history, r/mu grid value).
Fit: the r/mu minimising the unweighted sum of squared bin differences, with the simulated curves interpolated
linearly in log(r/mu) on a 2,000-point log grid between the end points. Uncertainty: 1,000 bootstrap resamples of
chromosomes (seed 20261030), simulated curves fixed; 95% percentile interval. Adequacy: at the fitted value the
simulated curve lies inside the observed bootstrap 95% band in >= 80% of bins. "shape" repeats the fit with every
curve divided by its first bin (sensitivity only).

Usage: 23e_fit_r_over_mu.py <obs_BAM.tsv> <obs_TP1.tsv> <ld_sim_dir> <out_fit.tsv> <out_curves.tsv>
"""
import csv, glob, os, sys
from collections import defaultdict
import numpy as np

NBIN = 50
NBOOT = 1000
SEED = 20261030
ADEQ = 0.80


def load_obs(fn):
    seg = defaultdict(lambda: [np.zeros(NBIN), np.zeros(NBIN)])
    for r in csv.DictReader(open(fn), delimiter="\t"):
        s = seg[(r["chrom"], r["block"])]
        s[0][int(r["bin"])] += int(r["n_pairs"])
        s[1][int(r["bin"])] += float(r["sum_r2"])
    by_chr = defaultdict(lambda: [np.zeros(NBIN), np.zeros(NBIN)])
    for (ch, _), (n, t) in seg.items():
        by_chr[ch][0] += n
        by_chr[ch][1] += t
    return by_chr


def load_sim(d):
    acc = defaultdict(lambda: [np.zeros(NBIN), np.zeros(NBIN)])
    for fn in glob.glob(os.path.join(d, "ld_*.tsv")):
        for r in csv.DictReader(open(fn), delimiter="\t"):
            a = acc[(r["model"], float(r["ratio"]), r["pop"])]
            a[0][int(r["bin"])] += int(r["n_pairs"])
            a[1][int(r["bin"])] += float(r["sum_r2"])
    curves = defaultdict(dict)
    for (m, k, pop), (n, t) in acc.items():
        curves[(m, pop)][k] = t / n
    return curves


def interp_family(grid):
    ks = np.array(sorted(grid))
    C = np.vstack([grid[k] for k in ks])
    lk = np.log(ks)
    fine = np.exp(np.linspace(lk[0], lk[-1], 2000))
    lf = np.log(fine)
    j = np.clip(np.searchsorted(lk, lf, "right") - 1, 0, len(ks) - 2)
    w = ((lf - lk[j]) / (lk[j + 1] - lk[j]))[:, None]
    return fine, (1 - w) * C[j] + w * C[j + 1]


def fit(obs, fine, F):
    sse = ((F - obs[None, :]) ** 2).sum(1)
    i = int(np.argmin(sse))
    return fine[i], i


def main():
    obs_fn = {"BAM": sys.argv[1], "TP1": sys.argv[2]}
    sim = load_sim(sys.argv[3])
    out_fit, out_curves = sys.argv[4], sys.argv[5]
    rng = np.random.default_rng(SEED)
    fit_rows, curve_rows = [], []
    for pop, fn in obs_fn.items():
        by_chr = load_obs(fn)
        chrs = sorted(by_chr)
        N = np.vstack([by_chr[c][0] for c in chrs])
        T = np.vstack([by_chr[c][1] for c in chrs])
        obs = T.sum(0) / N.sum(0)
        boots = []
        for _ in range(NBOOT):
            idx = rng.integers(0, len(chrs), len(chrs))
            boots.append(T[idx].sum(0) / N[idx].sum(0))
        boots = np.vstack(boots)
        lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
        for b in range(NBIN):
            curve_rows.append(dict(kind="observed", pop=pop, model="", r_over_mu="", bin=b,
                                   mid_kb=(b + 0.5) * 2, r2=obs[b], lo=lo[b], hi=hi[b]))
        for (model, spop), grid in sorted(sim.items()):
            if spop != pop:
                continue
            fine, F = interp_family(grid)
            for kind in ("full", "shape"):
                if kind == "full":
                    Fk, ok_, bk = F, obs, boots
                else:
                    Fk, ok_, bk = F / F[:, :1], obs / obs[0], boots / boots[:, :1]
                k_hat, i_hat = fit(ok_, fine, Fk)
                kb = np.array([fit(bk[j], fine, Fk)[0] for j in range(NBOOT)])
                klo, khi = np.percentile(kb, [2.5, 97.5])
                blo, bhi = np.percentile(bk, [2.5, 97.5], axis=0)
                inside = float(np.mean((Fk[i_hat] >= blo) & (Fk[i_hat] <= bhi)))
                fit_rows.append(dict(pop=pop, model=model, fit=kind, r_over_mu=k_hat, lo95=klo, hi95=khi,
                                     at_grid_edge=int(i_hat in (0, len(fine) - 1)),
                                     frac_bins_in_band=inside, adequate=int(inside >= ADEQ),
                                     sse=float(((Fk[i_hat] - ok_) ** 2).sum())))
                if kind == "full":
                    for k, c in sorted(grid.items()):
                        for b in range(NBIN):
                            curve_rows.append(dict(kind="simulated", pop=pop, model=model, r_over_mu=k, bin=b,
                                                   mid_kb=(b + 0.5) * 2, r2=c[b], lo="", hi=""))
    for fn, rows in ((out_fit, fit_rows), (out_curves, curve_rows)):
        with open(fn, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
    for r in fit_rows:
        print(f"{r['pop']} {r['model']:4s} {r['fit']:5s} r/mu {r['r_over_mu']:.3g} [{r['lo95']:.3g}, {r['hi95']:.3g}]"
              f" edge={r['at_grid_edge']} in-band {r['frac_bins_in_band']:.2f} adequate={r['adequate']}")


if __name__ == "__main__":
    main()
