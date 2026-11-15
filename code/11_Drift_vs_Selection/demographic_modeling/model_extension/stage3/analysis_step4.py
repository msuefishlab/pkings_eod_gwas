#!/usr/bin/env python3
"""analysis_step4.py -- false-rejection-rate (FRR) calibration of the criteria by simulation, one truth model.

Each of the 100 pseudo-observed M-sum genomes i (truth model, true r/mu) is pushed through the real-data pipeline:
  r/mu refit on its 13 calibration chromosomes against the truth model's calibration reference family (23e; BAM) ->
  r_i (+ edge flag); held-out band = its own 12-chromosome bootstrap band (ext_aggregate);
  P1 (held-out, gating form): simulated held-out curve at r_i vs i's held-out band, >= 80% of bins inside;
  P1 in-sample (continuity form): 23e on all 860 blocks against the in-sample family (genome band);
  P3: i's held-out BAM mean r^2 over 20-100 kb inside the central 95% of the replicate distribution at r_i;
  supporting diagnostics at 95% and 99%.
Implementation: 'replicates at r_i' are realised as
  - curves at r_i: the truth model's held-out scoring family (5 x 392 held-out block-sims per grid value) interpolated at
    r_i with the 23e log-linear rule;
  - replicate distributions at r_i: the other 99 pseudo-observed genomes (leave-one-out), for r/mu-dependent LD
    statistics shifted by m(r_i) - m(r_true), m = the scoring family's pooled statistic interpolated in log r/mu; SFS,
    divergence and interchromosomal statistics are not shifted;
  - Delta ladder: O/E against [1/3, 3] with E = mean of the other 99 genomes (unfolded truth models only).
Usage: analysis_step4.py <genomes_dir/<model>/r30> <fam_summary.npz> <chrom_split.tsv> <model> <r_true> <analysis>
                         <out_prefix>
"""
import csv, glob, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analysis_lib as A  # noqa: E402

gd, fsn, splitf, model, r_true, analysis, outp = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], \
    float(sys.argv[5]), sys.argv[6], sys.argv[7]
LEVEL = os.environ.get("S3_P1_LEVEL", "loo")    # 'loo' (adopted) or 'family' (first implementation, reported)
fs = dict(np.load(fsn)); split = A.read_split(splitf)
gdirs = sorted(glob.glob(os.path.join(gd, "g*")), key=lambda p: int(os.path.basename(p)[1:]))
G = [int(os.path.basename(p)[1:]) for p in gdirs]
n = len(gdirs)
fam_cal = {p: A.family(fs, model, "cal", "primary", p) for p in ("BAM", "TP1")}
fam_ins = {p: A.family(fs, model, "ins", "primary", p) for p in ("BAM", "TP1")}
fam_ho = {p: A.family(fs, model, "ho", "primary", p) for p in ("BAM", "TP1")}
far_ho = A.family_far(fs, model, "primary", "BAM")
b0_ho = A.family_far(fs, model, "primary", "BAM", bins=slice(0, 1))
cal_c, gen_c = A.set_chroms(split, "calibration"), A.set_chroms(split, "genome")

DIAG95 = [("pi", "BAM"), ("pi", "TP1"), ("dxy", "BAM|TP1"), ("fst_hudson", "BAM|TP1"), ("interchrom_r2_excess", "BAM"),
          ("interchrom_r2_excess", "TP1")]
SHIFTED = [("ld_r2_0_2kb", "BAM")]
DELTA = [(f"delta_ge_{t}_count", "BAM|TP1") for t in ("0.5", "0.7", "0.78", "0.8")]
per = []
panels = [A.read_panel(os.path.join(p, "panel.tsv")) for p in gdirs]
for i, p in enumerate(gdirs):
    z = os.path.join(p, "chromld.npz")
    cur = A.read_curves(os.path.join(p, "ld_curves.tsv"))
    rec = dict(g=G[i])
    for pop in ("BAM", "TP1"):
        N, T = A.chrom_sums(z, pop, cal_c)
        f = A.fit23e(N, T, fam_cal[pop], boots=(pop == "BAM"))
        rec[f"{pop}_r_hat"], rec[f"{pop}_edge"] = f["r_hat"], f["edge"]
        if pop == "BAM":
            rec["BAM_lo95"], rec["BAM_hi95"] = f["lo95"], f["hi95"]
        h = cur[(pop, "heldout")]
        rec[f"{pop}_halfband_heldout_mean"] = float(np.mean((h["hi"] - h["lo"]) / 2))
        if LEVEL == "family":
            fr, ab, be = A.bins_in_band(A.interp_at(fam_ho[pop], f["r_hat"]), h["lo"], h["hi"])
            rec[f"{pop}_P1_heldout_frac"], rec[f"{pop}_P1_heldout_above"], rec[f"{pop}_P1_heldout_below"] = fr, ab, be
            rec[f"{pop}_P1_heldout_pass"] = int(fr >= 0.8)
        N, T = A.chrom_sums(z, pop, gen_c)
        fi = A.fit23e(N, T, fam_ins[pop], boots=False)
        rec[f"{pop}_P1_insample_r_hat"], rec[f"{pop}_P1_insample_frac"] = fi["r_hat"], fi["frac_in_band_own"]
        rec[f"{pop}_P1_insample_pass"] = int(fi["frac_in_band_own"] >= 0.8)
        if pop == "BAM":
            fsh = A.fit23e(N, T, fam_ins[pop], shape=True, boots=False)
            rec["BAM_shape_insample_frac"] = fsh["frac_in_band_own"]
            rec["BAM_shape_insample_pass"] = int(fsh["frac_in_band_own"] >= 0.8)
    per.append(rec)
# held-out pooled sums per genome (for leave-one-out fresh-replicate curves)
ho_c = A.set_chroms(split, "heldout")
HS = {pop: np.array([np.vstack(A.chrom_sums(os.path.join(p, "chromld.npz"), pop, ho_c)).reshape(2, len(ho_c), 50).sum(1)
                     for p in gdirs]) for pop in ("BAM", "TP1")}          # (n, 2 [N, T], 50)
lk = np.log(A.GRID)


def slope(fam, r):
    """local d(curve)/d(log r/mu) of the scoring family on the grid interval containing r (23e log-linear rule)"""
    j = int(np.clip(np.searchsorted(lk, np.log(r), "right") - 1, 0, len(lk) - 2))
    return (np.asarray(fam[A.GRID[j + 1]]) - np.asarray(fam[A.GRID[j]])) / (lk[j + 1] - lk[j])


# P1 (LOO level + local slope), P3 and supporting diagnostics, leave-one-out
for i, rec in enumerate(per):
    others = [j for j in range(n) if j != i]
    if LEVEL == "loo":
        cur = A.read_curves(os.path.join(gdirs[i], "ld_curves.tsv"))
        for pop in ("BAM", "TP1"):
            Nn = HS[pop][others, 0].sum(0); Tt = HS[pop][others, 1].sum(0)
            curve = Tt / Nn + slope(fam_ho[pop], r_true) * (np.log(rec[f"{pop}_r_hat"]) - np.log(r_true))
            h = cur[(pop, "heldout")]
            fr, ab, be = A.bins_in_band(curve, h["lo"], h["hi"])
            rec[f"{pop}_P1_heldout_frac"], rec[f"{pop}_P1_heldout_above"], rec[f"{pop}_P1_heldout_below"] = fr, ab, be
            rec[f"{pop}_P1_heldout_pass"] = int(fr >= 0.8)
    dl = np.log(rec["BAM_r_hat"]) - np.log(r_true)
    d_far = float(slope(far_ho, r_true) * dl) if LEVEL == "loo" else A.interp_at(far_ho, rec["BAM_r_hat"]) - A.interp_at(far_ho, r_true)
    d_b0 = float(slope(b0_ho, r_true) * dl) if LEVEL == "loo" else A.interp_at(b0_ho, rec["BAM_r_hat"]) - A.interp_at(b0_ho, r_true)
    x = panels[i][("ld_r2_mean_20_100kb", "BAM", "heldout")][0]
    ref = np.array([panels[j][("ld_r2_mean_20_100kb", "BAM", "heldout")][0] for j in others]) + d_far
    lo, hi = np.percentile(ref, [2.5, 97.5])
    rec["P3_value"], rec["P3_lo"], rec["P3_hi"], rec["P3_pass"] = x, lo, hi, int(lo <= x <= hi)
    rec["P1_and_P3_pass"] = int(rec["BAM_P1_heldout_pass"] and rec["P3_pass"])
    for stat, smp in DIAG95 + SHIFTED:
        for s in ("genome", "heldout"):
            x = panels[i][(stat, smp, s)][0]
            ref = np.array([panels[j][(stat, smp, s)][0] for j in others]) + (d_b0 if (stat, smp) in SHIFTED else 0)
            for lev, q in ((95, (2.5, 97.5)), (99, (0.5, 99.5))):
                lo, hi = np.percentile(ref, q)
                rec[f"{stat}|{smp}|{s}|pass{lev}"] = int(lo <= x <= hi)
    if analysis == "unfolded":
        for stat, smp in DELTA:
            for s in ("genome", "heldout"):
                x = panels[i][(stat, smp, s)][0]
                E = np.mean([panels[j][(stat, smp, s)][0] for j in others])
                oe = x / E if E > 0 else np.nan
                rec[f"{stat}|{smp}|{s}|OE"] = oe
                rec[f"{stat}|{smp}|{s}|passOE"] = int(1 / 3 <= oe <= 3)
cols = list(per[0])
with open(outp + "_pseudo_obs.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(per)
# FRR table
frr = []
crit = [("P1_heldout (BAM, gating form)", "BAM_P1_heldout_pass"), ("P3 (BAM held-out 20-100 kb)", "P3_pass"),
        ("P1_heldout AND P3", "P1_and_P3_pass"), ("P1_insample (BAM, continuity form)", "BAM_P1_insample_pass"),
        ("TP1 curve, P1 rule (held-out)", "TP1_P1_heldout_pass"), ("TP1 curve, in-sample", "TP1_P1_insample_pass"),
        ("BAM shape-only fit (in-sample)", "BAM_shape_insample_pass")]
crit += [(k, k) for k in cols if k.endswith("|pass95") or k.endswith("|pass99") or k.endswith("|passOE")]
for lab, k in crit:
    fails = int(sum(1 - r[k] for r in per))
    lo, hi = A.clopper_pearson(fails, n)
    frr.append(dict(truth_model=model, r_true=r_true, criterion=lab, n=n, n_fail=fails, FRR=fails / n, CP_lo=lo, CP_hi=hi))
with open(outp + "_frr.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=list(frr[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(frr)
rh = np.array([r["BAM_r_hat"] for r in per])
print(f"{model} r_true={r_true} n={n} refit r/mu median {np.median(rh):.3g} [{np.percentile(rh, 2.5):.3g}, "
      f"{np.percentile(rh, 97.5):.3g}] edge {int(sum(r['BAM_edge'] for r in per))}")
for r in frr[:7]:
    print(f"  {r['criterion']:40s} FRR {r['FRR']:.2f} [{r['CP_lo']:.3f}, {r['CP_hi']:.3f}]")
