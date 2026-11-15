#!/usr/bin/env python3
"""ext_check_23a.py -- does the panel's LD reproduce simulation/23a per-block rows and the 23e observed
curve?

1. Block level: for the 200 blocks of 23a, panel ld_cnt / ld_tot / n_ld_snps vs 23a_observed_ld_<pop>.tsv.
2. Curve level: 23e's observed curve and band recomputed from the panel's values for those 200 blocks with 23e's
   own bootstrap (one rng(20261030), BAM then TP1, chromosomes sorted) vs 23_ld_curves.tsv observed rows.
3. Full mask (all eligible blocks) vs 23_ld_curves.tsv observed rows: difference per bin, share inside 23e band.
Usage: ext_check_23a.py <seg_dir> <simulation_output_dir> <out.tsv>
"""
import csv, os, sys
from collections import defaultdict
import numpy as np

seg_dir, sim, out = sys.argv[1:4]
NB = 50
blk = defaultdict(lambda: np.zeros(NB))
sc = {}
for fn in os.listdir(seg_dir):
    if fn.startswith("seg.") and fn.endswith(".tsv"):
        for line in open(os.path.join(seg_dir, fn)):
            ch, b, kind, name, i, v = line.rstrip("\n").split("\t")
            if name in ("BAM", "TP1") and kind in ("ld_cnt", "ld_tot"):
                blk[(name, kind, ch, int(b))][int(i)] = float(v)
            elif name in ("BAM", "TP1") and kind == "n_ld_snps":
                sc[(name, ch, int(b))] = int(float(v))
rows = []
curves = {(r["pop"], int(r["bin"])): r for r in csv.DictReader(open(os.path.join(sim, "23_ld_curves.tsv")), delimiter="\t")
          if r["kind"] == "observed"}
rng = np.random.default_rng(20261030)
for pop in ("BAM", "TP1"):
    ref = defaultdict(lambda: [np.zeros(NB), np.zeros(NB), 0])
    for r in csv.DictReader(open(os.path.join(sim, f"23a_observed_ld_{pop}.tsv")), delimiter="\t"):
        k = (r["chrom"], int(r["block"]))
        ref[k][0][int(r["bin"])] = int(r["n_pairs"]); ref[k][1][int(r["bin"])] = float(r["sum_r2"]); ref[k][2] = int(r["n_snps"])
    dn = max(np.abs(blk[(pop, "ld_cnt") + k] - v[0]).max() for k, v in ref.items())
    dt = max(np.abs(blk[(pop, "ld_tot") + k] - v[1]).max() for k, v in ref.items())
    ds = sum(sc.get((pop,) + k, -1) != v[2] for k, v in ref.items())
    rows.append(dict(check="block_max_abs_diff_n_pairs", pop=pop, value=dn))
    rows.append(dict(check="block_max_abs_diff_sum_r2", pop=pop, value=dt))
    rows.append(dict(check="blocks_with_n_snps_mismatch_of_200", pop=pop, value=ds))
    byc = defaultdict(lambda: [np.zeros(NB), np.zeros(NB)])
    for (ch, b) in ref:
        byc[ch][0] += blk[(pop, "ld_cnt", ch, b)]; byc[ch][1] += blk[(pop, "ld_tot", ch, b)]
    chrs = sorted(byc)
    N = np.vstack([byc[c][0] for c in chrs]); T = np.vstack([byc[c][1] for c in chrs])
    obs = T.sum(0) / N.sum(0)
    boots = []
    for _ in range(1000):
        idx = rng.integers(0, len(chrs), len(chrs))
        boots.append(T[idx].sum(0) / N[idx].sum(0))
    lo, hi = np.percentile(np.vstack(boots), [2.5, 97.5], axis=0)
    r0 = np.array([float(curves[(pop, b)]["r2"]) for b in range(NB)])
    l0 = np.array([float(curves[(pop, b)]["lo"]) for b in range(NB)])
    h0 = np.array([float(curves[(pop, b)]["hi"]) for b in range(NB)])
    rows.append(dict(check="curve200_max_abs_diff_r2", pop=pop, value=np.abs(obs - r0).max()))
    rows.append(dict(check="curve200_max_abs_diff_lo", pop=pop, value=np.abs(lo - l0).max()))
    rows.append(dict(check="curve200_max_abs_diff_hi", pop=pop, value=np.abs(hi - h0).max()))
    allc = defaultdict(lambda: [np.zeros(NB), np.zeros(NB)])
    for (p_, kind, ch, b), v in blk.items():
        if p_ == pop:
            allc[kind][0] += v
    full = allc["ld_tot"][0] / allc["ld_cnt"][0]
    rows.append(dict(check="fullmask_max_abs_diff_r2_vs_23e", pop=pop, value=np.abs(full - r0).max()))
    rows.append(dict(check="fullmask_mean_diff_r2_vs_23e", pop=pop, value=(full - r0).mean()))
    rows.append(dict(check="fullmask_share_bins_inside_23e_band", pop=pop, value=np.mean((full >= l0) & (full <= h0))))
    rows.append(dict(check="fullmask_n_blocks", pop=pop, value=len({k[2:] for k in blk if k[0] == pop})))
with open(out, "w") as fh:
    w = csv.DictWriter(fh, fieldnames=["check", "pop", "value"], delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
for r in rows:
    print(r)
