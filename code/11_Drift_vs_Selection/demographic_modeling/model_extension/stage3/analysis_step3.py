#!/usr/bin/env python3
"""analysis_step3.py -- criterion P1 r/mu fits per model (physical r/mu, the 12-value grid of simulation/23e).
  calibration form : observed calibration-chromosome curve (13 chr) vs the model's 'cal' reference family -> fitted
                     r/mu, 95% bootstrap interval (its lower bound = the 23e lower bound), grid-edge flag; BAM under
                     every treatment row (callable / d019 against their own observed rows; gqdrop / eps against the
                     primary observed curve); TP1 primary (supporting diagnostic).
  in-sample form   : simulation/23e as run (all 860 M-sum blocks, 'ins' family, interpolated curve at the fit vs the
                     genome band, >= 80% of bins) -- continuity; 'shape' = 23e shape-only fit (supporting).
Usage: analysis_step3.py <fam_summary.npz> <obs_dir (obs_<row>_chromld.npz, obs_<row>_ld_curves.tsv)> <chrom_split.tsv>
                         <models,comma> <out_fits.tsv>
"""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analysis_lib as A  # noqa: E402

fsn, od, splitf, models, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(","), sys.argv[5]
fs = dict(np.load(fsn))
split = A.read_split(splitf)
OBSROW = {"primary": "primary", "callable": "callable", "gqdrop": "primary", "d019": "d019", "eps0.005": "primary",
          "eps0.01": "primary"}
rows, chk = [], []
cache = {}


def obs(row, pop, s):
    k = (row, pop, s)
    if k not in cache:
        cs = A.set_chroms(split, s)
        cache[k] = A.chrom_sums(os.path.join(od, f"obs_{row}_chromld.npz"), pop, cs)
    return cache[k]


# check: bands recomputed here == ext_aggregate bands of the observed rows
for row in ("primary", "callable", "d019"):
    cur = A.read_curves(os.path.join(od, f"obs_{row}_ld_curves.tsv"))
    for pop in ("BAM", "TP1"):
        for s in ("genome", "calibration", "heldout"):
            N, T = obs(row, pop, s)
            o, bt, lo, hi = A.obs_curve_boot(N, T)
            c = cur[(pop, s)]
            chk.append((row, pop, s, float(np.max(np.abs(o - c["r2"]))), float(np.max(np.abs(lo - c["lo"]))),
                        float(np.max(np.abs(hi - c["hi"])))))
for m in models:
    for t, row in OBSROW.items():
        N, T = obs(row, "BAM", "calibration")
        r = A.fit23e(N, T, A.family(fs, m, "cal", t, "BAM"))
        rows.append(dict(model=m, pop="BAM", treatment=t, obs_row=row, form="calibration", **{k: r[k] for k in (
            "r_hat", "lo95", "hi95", "edge", "frac_in_band_own", "sse")}))
    N, T = obs("primary", "TP1", "calibration")
    r = A.fit23e(N, T, A.family(fs, m, "cal", "primary", "TP1"))
    rows.append(dict(model=m, pop="TP1", treatment="primary", obs_row="primary", form="calibration",
                     **{k: r[k] for k in ("r_hat", "lo95", "hi95", "edge", "frac_in_band_own", "sse")}))
    for pop in ("BAM", "TP1"):
        N, T = obs("primary", pop, "genome")
        for form, shape in (("in_sample", False), ("in_sample_shape", True)):
            r = A.fit23e(N, T, A.family(fs, m, "ins", "primary", pop), shape=shape)
            rows.append(dict(model=m, pop=pop, treatment="primary", obs_row="primary", form=form,
                             **{k: r[k] for k in ("r_hat", "lo95", "hi95", "edge", "frac_in_band_own", "sse")}))
for r in rows:
    r["adequate_own_80"] = int(r["frac_in_band_own"] >= 0.8)
with open(out + ".tmp", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
os.replace(out + ".tmp", out)
with open(out.replace(".tsv", "_bandcheck.tsv"), "w") as fo:
    fo.write("row\tpop\tset\tmax_abs_d_r2\tmax_abs_d_lo\tmax_abs_d_hi\n")
    for x in chk:
        fo.write("\t".join(str(v) for v in x) + "\n")
print("fits", len(rows), "max band discrepancy", max(max(x[3:]) for x in chk))
