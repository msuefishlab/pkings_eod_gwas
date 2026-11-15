#!/usr/bin/env python3
"""pilot_compare.py -- Stage 3 pipeline validation: this pipeline's re-simulation of the 40 coord-pilot blocks under M2
(pilot seeds) vs the pilot's arm b'' task files (coord_pilot/a2/sim/bpp_<ib>_<ir>.npz), per block x r/mu x population x
scoring (bpp = primary LD site set, b_chk = callable subset): exact equality of n_pairs and sum r^2 per replicate and
bin, and the pooled-curve difference in units of the Stage 1 held-out half-band, with the MC SE of the pooled pilot
curve (bootstrap over replicate index, 200 resamples, rng 20261112) for the 'within MC error' statement.
Usage: pilot_compare.py <new_dir> <pilot_sim_dir> <observed_ld_curves.tsv> <out.tsv>
"""
import csv, os, sys
import numpy as np

new, old, obs, out = sys.argv[1:5]
hb = {}
for r in csv.DictReader(open(obs), delimiter="\t"):
    if r["kind"] == "observed" and r["set"] == "heldout" and r["pop"] in ("BAM", "TP1"):
        hb[(r["pop"], int(r["bin"]))] = (float(r["hi"]) - float(r["lo"])) / 2
rows, allc = [], {}
for ir in range(5):
    for ib in range(40):
        a = np.load(os.path.join(new, f"pilot_{ib}_{ir}.npz")); b = np.load(os.path.join(old, f"bpp_{ib}_{ir}.npz"))
        for v in ("bpp", "b_chk"):
            for pop in ("BAM", "TP1"):
                cn, tn = a[f"{v}|ld_cnt|{pop}"], a[f"{v}|ld_tot|{pop}"]
                co, to = b[f"{v}|ld_cnt|{pop}"], b[f"{v}|ld_tot|{pop}"]
                eq = bool(np.array_equal(cn, co) and np.array_equal(tn, to))
                rows.append(dict(i_ratio=ir, i_block=ib, scoring=v, pop=pop, identical=int(eq),
                                 max_abs_dtot=float(np.max(np.abs(tn - to))) if tn.shape == to.shape else float("nan")))
                k = (ir, v, pop)
                allc.setdefault(k, [[], [], [], []])
                for L_, x in zip(allc[k], (cn, tn, co, to)):
                    L_.append(x)
summ = []
rng = np.random.default_rng(20261112)
for (ir, v, pop), (cn, tn, co, to) in sorted(allc.items()):
    CN, TN, CO, TO = (np.sum(x, axis=0) for x in (cn, tn, co, to))        # (100 reps, 50) summed over blocks
    new_c = TN.sum(0) / CN.sum(0); old_c = TO.sum(0) / CO.sum(0)
    bs = [];
    for _ in range(200):
        i = rng.integers(0, CO.shape[0], CO.shape[0]); bs.append(TO[i].sum(0) / CO[i].sum(0))
    se = np.std(bs, axis=0)
    h = np.array([hb[(pop, b)] for b in range(50)])
    summ.append(dict(i_ratio=ir, scoring=v, pop=pop, max_abs_diff_halfbands=float(np.max(np.abs(new_c - old_c) / h)),
                     max_mc_se_halfbands=float(np.max(se / h)),
                     within_2se=int(np.all(np.abs(new_c - old_c) <= 2 * se + 1e-15))))
with open(out + ".tmp", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(summ[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(summ)
os.replace(out + ".tmp", out)
n_id = sum(r["identical"] for r in rows)
print(f"identical {n_id}/{len(rows)} (block x r/mu x scoring x pop); all within 2 MC SE:",
      all(s["within_2se"] for s in summ), flush=True)
with open(out.replace(".tsv", "_detail.tsv"), "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
