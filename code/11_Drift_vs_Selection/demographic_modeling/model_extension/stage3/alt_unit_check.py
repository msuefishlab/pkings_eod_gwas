#!/usr/bin/env python3
"""alt_unit_check.py -- unit check of the alternative-model conversion: simulated joint SFS at callable
positions (drift-test projection 60 x 60; sim_task.py sfs) vs dadi's expected spectrum at the fit (alt_expected_sfs.py),
scaled by L_sim / L_mask. Folding and masking use dadi.Spectrum (fold() for folded models; the data mask).
Reported per model: simulated / expected for (a) all unmasked cells, (b) each population's marginal spectrum in 6
frequency classes, (c) the Delta-ladder classes (D >= 0.5, 0.7, 0.8 unfolded; |D| folded); z = (sim - exp) / SE,
SE = bootstrap over simulated units (blocks; 500 resamples, rng 20261114). Run inside dadi_tools.sif.
Usage: alt_unit_check.py <expected.npz> <sfs_dir> <models,comma> <out.tsv>
"""
import csv, glob, os, sys
import numpy as np
import dadi

expf, sd, models, out = sys.argv[1], sys.argv[2], sys.argv[3].split(","), sys.argv[4]
E = np.load(expf)
rows = []
rng = np.random.default_rng(20261114)
n1 = n2 = 60
I, J = np.meshgrid(np.arange(n1 + 1), np.arange(n2 + 1), indexing="ij")
D = I / n1 - J / n2
CL = [(1, 2), (3, 6), (7, 15), (16, 30), (31, 45), (46, 59)]
for m in models:
    exp, mask, Lm = E[f"{m}|exp"], E[f"{m}|mask"], float(E[f"{m}|L_mask"])
    folded = m.endswith("f")
    fl = sorted(glob.glob(os.path.join(sd, f"{m}_*.npz")))
    per = np.concatenate([np.load(f)["per_fs"] for f in fl]); pl = np.concatenate([np.load(f)["per_L"] for f in fl])

    def prep(F):
        s = dadi.Spectrum(F)
        if folded:
            s = s.fold()
        return np.ma.filled(np.ma.masked_array(np.ma.filled(s, 0.0), mask), 0.0)
    units = np.array([prep(F) for F in per])
    Ls = pl.sum()
    ex = np.where(mask, 0.0, exp) * Ls / Lm
    stats = {"all unmasked cells": np.ones_like(D, bool)}
    for a, b in CL:
        stats[f"BAM marginal {a}-{b}/60"] = (I >= a) & (I <= b)
        stats[f"TP1 marginal {a}-{b}/60"] = (J >= a) & (J <= b)
    for t in (0.5, 0.7, 0.8):
        stats[f"{'|D|' if folded else 'D'} >= {t}"] = (np.abs(D) if folded else D) >= t - 1e-9
    for lab, sel in stats.items():
        s_units = units[:, sel].sum(1)
        sim = s_units.sum(); e = ex[sel].sum()
        bs = [s_units[rng.integers(0, len(s_units), len(s_units))].sum() for _ in range(500)]
        se = float(np.std(bs))
        rows.append(dict(model=m, statistic=lab, simulated=sim, expected=e, ratio=sim / e if e > 0 else float("nan"),
                         se=se, z=(sim - e) / se if se > 0 else float("nan"), n_units=len(units), L_sim=int(Ls)))
with open(out, "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
for m in models:
    rr = [r for r in rows if r["model"] == m]
    print(m, "all-cells ratio %.4f z %.2f; max |z| %.2f" % (rr[0]["ratio"], rr[0]["z"], max(abs(r["z"]) for r in rr)))
