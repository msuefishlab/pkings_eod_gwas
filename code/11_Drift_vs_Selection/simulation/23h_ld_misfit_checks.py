#!/usr/bin/env python3
"""23h_ld_misfit_checks.py -- why does no fitted history reproduce the observed LD decay? (LD-decay diagnostic)

Three checks on existing outputs (23a, 23e), no new data:
  1. Mixture: can a mix of the simulated uniform-rate curves of one history (weights >= 0, summing to 1) reproduce
     the observed mean curve? Tests recombination-rate variation, within the r/mu grid.
  2. Spread: long-range r^2 (80-100 kb, layout units) per 1-Mb segment, against the simulated mean (M2). If a few
     low-recombination segments drove the misfit, most segments would sit near the simulation.
  3. Chromosome position: does long-range r^2 rise away from chromosome ends? Expected if crossovers concentrate
     near the ends, as in many teleosts.

Usage (repo root): 23h_ld_misfit_checks.py <fai> > 23h_ld_misfit_checks.txt
"""
import csv, sys
from collections import defaultdict
import numpy as np
from scipy.optimize import nnls
from scipy.stats import spearmanr

S = "output_data/11_Drift_vs_Selection/simulation"
FAR = slice(40, 50)                    # 80-100 kb
L = {l.split()[0]: int(l.split()[1]) for l in open(sys.argv[1])}
rows = list(csv.DictReader(open(f"{S}/23_ld_curves.tsv"), delimiter="\t"))


def sim_curves(pop, model):
    cur = {}
    for r in rows:
        if r["kind"] == "simulated" and r["pop"] == pop and r["model"] == model:
            cur.setdefault(float(r["r_over_mu"]), np.zeros(50))[int(r["bin"])] = float(r["r2"])
    return cur


for pop in ("BAM", "TP1"):
    o = sorted([r for r in rows if r["kind"] == "observed" and r["pop"] == pop], key=lambda r: int(r["bin"]))
    obs, lo, hi = (np.array([float(r[c]) for r in o]) for c in ("r2", "lo", "hi"))
    print(f"== {pop}")
    for m in sorted({r["model"] for r in rows if r["kind"] == "simulated"}):
        cur = sim_curves(pop, m)
        ks = sorted(cur)
        C = np.vstack([cur[k] for k in ks]).T
        w, _ = nnls(np.vstack([C, 1e3 * np.ones(len(ks))]), np.append(obs, 1e3))
        fit = C @ w
        mix = ", ".join(f"{k:g}:{x:.2f}" for k, x in zip(ks, w) if x > 0.005)
        print(f"  mixture {m:3s} bins in band {np.mean((fit >= lo) & (fit <= hi)):.2f}  weights {mix}")
    seg = defaultdict(lambda: [np.zeros(50), np.zeros(50)])
    for r in csv.DictReader(open(f"{S}/23a_observed_ld_{pop}.tsv"), delimiter="\t"):
        s = seg[(r["chrom"], int(r["block"]))]
        s[0][int(r["bin"])] += int(r["n_pairs"])
        s[1][int(r["bin"])] += float(r["sum_r2"])
    x, y = [], []
    for (ch, b), (n, t) in seg.items():
        if n[FAR].sum() == 0:
            continue
        mid = (b + 0.5) * 1e6
        x.append(min(mid, L[ch] - mid) / L[ch])
        y.append(t[FAR].sum() / n[FAR].sum())
    x, y = np.array(x), np.array(y)
    cur = sim_curves(pop, "M2")
    s025, s1 = cur[0.25][FAR].mean(), cur[1.0][FAR].mean()
    print(f"  spread: {y.size} segments, r2 at 80-100 kb p10 {np.percentile(y, 10):.3f} median {np.median(y):.3f}"
          f" p90 {np.percentile(y, 90):.3f}; M2 simulated mean r/mu 0.25 {s025:.3f}, r/mu 1 {s1:.3f}")
    print(f"  share of segments above M2 r/mu 1: {np.mean(y > s1):.2f}; above M2 r/mu 0.25: {np.mean(y > s025):.2f}")
    rho, p = spearmanr(x, y)
    print(f"  position: Spearman(distance from nearest end / length, r2 at 80-100 kb) = {rho:.2f}, p = {p:.2g}")
    for a, b in ((0, 0.1), (0.1, 0.25), (0.25, 0.5)):
        k = (x >= a) & (x < b)
        print(f"    {a:.2f}-{b:.2f} from end: n {k.sum():3d}, median r2 {np.median(y[k]):.3f}")
