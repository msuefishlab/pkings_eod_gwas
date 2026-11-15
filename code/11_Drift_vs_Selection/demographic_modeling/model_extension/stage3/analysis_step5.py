#!/usr/bin/env python3
"""analysis_step5.py -- scoring of the real data for one model (run only after the FRR table of analysis_step4.py
is fixed). Replicates: the model's 100 fresh whole M-sum genomes at its BAM calibration-fitted r/mu (i_ratio 20).
  P1 held-out (gating): pooled simulated held-out BAM curve (100 genomes x 392 blocks) vs the observed held-out band
     (Stage 1, 12-chromosome bootstrap); pass iff >= 80% of bins inside. Grid-edge fits are reported as bounds.
  P1 in-sample (continuity): simulation/23e as run (from analysis_step3 fits).
  P2: drift-test admissibility of the dadi fit (drift-test / Stage 2b table); not re-simulated.
  P3: observed held-out BAM mean r^2 20-100 kb inside the central 95% of the 100 genomes' values; gating iff the
     combined FRR(P1 and P3) <= 0.10 under every truth model scored, else a supporting diagnostic.
  Supporting: central 95% (99% where fallback ii applies) of the replicate distribution, genome-wide and held-out;
     Delta ladder as O/E in [1/3, 3] (unfolded models); TP1 curve under the P1 rule (TP1 calibration fit, held-out
     scoring family at that r/mu); BAM shape-only fit (in-sample 23e); recent BAM Ne (last 10 generations, both mu arms)
     against the temporal-Ne lower bounds (>= ~3,100 at 1 yr/generation, >= ~1,500 at 2 yr).
  Sensitivity rows (never gating): fitted r/mu (calibration) and held-out bins-in-band for callable / gqdrop / d019 /
     eps0.005 / eps0.01, bins-in-band from the held-out scoring family at the row's fitted r/mu; shifts are relative to
     the primary row evaluated the same way. Fixed points: primary held-out bins-in-band at physical r/mu
     0.25, 1.0, 2.5 from the scoring family (exact grid values).
Usage: analysis_step5.py <model> <genomes_dir/<model>/r20> <fam_summary.npz> <fits.tsv> <obs3_dir> <stage1_obs_panel>
         <chrom_split> <params.json> <analysis> <admissible 0|1> <p3_gating 0|1> <diag99 list or -> <out_prefix>
"""
import csv, glob, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analysis_lib as A  # noqa: E402

(model, gd, fsn, fitsf, od, s1p, splitf, pj, analysis, admiss, p3g, d99, outp) = sys.argv[1:14]
admiss, p3g = int(admiss), int(p3g)
D99 = set() if d99 == "-" else set(d99.split(","))
fs = dict(np.load(fsn)); split = A.read_split(splitf)
fits = [r for r in csv.DictReader(open(fitsf), delimiter="\t") if r["model"] == model]
fit = {(r["pop"], r["treatment"], r["form"]): r for r in fits}
gdirs = sorted(glob.glob(os.path.join(gd, "g*")), key=lambda p: int(os.path.basename(p)[1:]))
n = len(gdirs)
panels = [A.read_panel(os.path.join(p, "panel.tsv")) for p in gdirs]
ho_c, cal_c, gen_c = A.set_chroms(split, "heldout"), A.set_chroms(split, "calibration"), A.set_chroms(split, "genome")
obs = {row: A.read_curves(os.path.join(od, f"obs_{row}_ld_curves.tsv")) for row in ("primary", "callable", "d019")}
opan = A.read_panel(os.path.join(od, "obs_primary_panel.tsv"))
rows = []
COMP = ("composition limit: the primary row matches the observed LD SNP set by position, not by allele-level "
        "callability (94.7% of simulated vs 64.5% of observed ld23 SNPs callable)")


def add(criterion, sset, value, lo, hi, outcome, gating, note="", **kw):
    rows.append(dict(model=model, analysis=analysis, criterion=criterion, set=sset, value=value, lo=lo, hi=hi,
                     outcome=outcome, gating=gating, note=note, **kw))


def pooled(pop, chroms):
    N = sum(A.chrom_sums(os.path.join(p, "chromld.npz"), pop, chroms)[0].sum(0) for p in gdirs)
    T = sum(A.chrom_sums(os.path.join(p, "chromld.npz"), pop, chroms)[1].sum(0) for p in gdirs)
    return T / N


def misfit(curve, band):
    h = (band["hi"] - band["lo"]) / 2
    d = (curve - band["r2"]) / h
    return float(np.mean(d[A.FAR])), float(d[0])


f_cal = fit[("BAM", "primary", "calibration")]
r_hat = float(f_cal["r_hat"])
sim_ho = pooled("BAM", ho_c)
band = obs["primary"][("BAM", "heldout")]
fr, ab, be = A.bins_in_band(sim_ho, band["lo"], band["hi"])
mf_far, mf_0 = misfit(sim_ho, band)
p1 = int(fr >= 0.8)
add("P1 held-out (gating)", "heldout", fr, "", "", "pass" if p1 else "fail", "gating",
    f"r/mu(cal) {r_hat:.4g} [{float(f_cal['lo95']):.4g}, {float(f_cal['hi95']):.4g}]"
    f"{' AT GRID EDGE (bound)' if f_cal['edge'] == '1' else ''}; bins above/below band {ab}/{be}; "
    f"mean misfit 20-100 kb {mf_far:+.2f} half-bands, 0-2 kb {mf_0:+.2f}; {COMP}",
    r_over_mu=r_hat, r_lo95=f_cal["lo95"], r_hi95=f_cal["hi95"], edge=f_cal["edge"])
fi = fit[("BAM", "primary", "in_sample")]
p1i = int(fi["adequate_own_80"])
add("P1 in-sample (continuity)", "genome", float(fi["frac_in_band_own"]), "", "", "pass" if p1i else "fail",
    "reported", f"r/mu(in-sample) {float(fi['r_hat']):.4g} [{float(fi['lo95']):.4g}, {float(fi['hi95']):.4g}]"
    f"{' AT GRID EDGE' if fi['edge'] == '1' else ''}"
    + ("" if p1i == p1 else "; LD adequacy depends on calibration/validation split"), r_over_mu=fi["r_hat"],
    r_lo95=fi["lo95"], r_hi95=fi["hi95"], edge=fi["edge"])
add("P2 SFS admissibility (drift-test rule)", "fit", admiss, "", "", "pass" if admiss else "fail", "gating", "dadi fit property")
x = opan[("ld_r2_mean_20_100kb", "BAM", "heldout")][0]
ref = np.array([pp[("ld_r2_mean_20_100kb", "BAM", "heldout")][0] for pp in panels])
lo, hi = np.percentile(ref, [2.5, 97.5])
p3 = int(lo <= x <= hi)
add("P3 long-range LD level", "heldout", x, lo, hi, ("pass" if p3 else ("fail: observed above" if x > hi else
                                                                       "fail: observed below")),
    "gating" if p3g else "supporting (demoted by the FRR calibration)", f"simulated median {np.median(ref):.4g}; {COMP}")
# supporting diagnostics
DIAG = [("pi", "BAM"), ("pi", "TP1"), ("dxy", "BAM|TP1"), ("fst_hudson", "BAM|TP1"), ("interchrom_r2_excess", "BAM"),
        ("ld_r2_0_2kb", "BAM"), ("interchrom_r2_excess", "TP1")]
s1 = {}
for r in csv.DictReader(open(s1p), delimiter="\t"):
    s1.setdefault((r["statistic"], r["sample"], r["set"]), float(r["value"]))
misfits = [] if (p3 or p3g) else ["P3 long-range LD level (demoted)"]
for stat, smp in DIAG:
    for s in ("genome", "heldout"):
        key = f"{stat}|{smp}|{s}"
        lev = 99 if key in D99 else 95
        q = (0.5, 99.5) if lev == 99 else (2.5, 97.5)
        x = opan[(stat, smp, s)][0]
        ref = np.array([pp[(stat, smp, s)][0] for pp in panels])
        lo, hi = np.percentile(ref, q)
        ok = lo <= x <= hi
        note = f"central {lev}%" + (" (fallback ii)" if lev == 99 else "")
        if stat == "interchrom_r2_excess" and smp == "TP1":
            note += ("; structurally unmatched (TP1 pools two localities); observed TP_APA2009 "
                     f"{s1.get((stat, 'TP_APA2009', s), float('nan')):.3g}, TP_BEN {s1.get((stat, 'TP_BEN', s), float('nan')):.3g}")
        out = "inside" if ok else ("outside: observed above" if x > hi else "outside: observed below")
        add(f"{stat} {smp}", s, x, lo, hi, out, "supporting", note)
        if not ok and not (stat == "interchrom_r2_excess" and smp == "TP1"):
            misfits.append(f"{stat} {smp} ({s})")
if analysis == "unfolded":
    for t in ("0.5", "0.7", "0.78", "0.8"):
        for s in ("genome", "heldout"):
            x = opan[(f"delta_ge_{t}_count", "BAM|TP1", s)][0]
            E = float(np.mean([pp[(f"delta_ge_{t}_count", "BAM|TP1", s)][0] for pp in panels]))
            oe = x / E if E > 0 else float("nan")
            ok = 1 / 3 <= oe <= 3
            add(f"Delta >= {t} O/E", s, oe, 1 / 3, 3, "inside" if ok else "outside", "supporting",
                f"O {x:g}, E {E:.4g} (mean of {n} genomes)")
            if not ok:
                misfits.append(f"Delta>={t} O/E ({s})")
else:
    add("Delta ladder O/E", "-", "", "", "", "not applicable (folded model: polarized statistic)", "supporting")
# TP1 curve under the P1 rule
ft = fit[("TP1", "primary", "calibration")]
cT = A.interp_at(A.family(fs, model, "ho", "primary", "TP1"), float(ft["r_hat"]))
bT = obs["primary"][("TP1", "heldout")]
frT, abT, beT = A.bins_in_band(cT, bT["lo"], bT["hi"])
mT = misfit(cT, bT)
add("TP1 LD curve, P1 rule", "heldout", frT, "", "", "inside (>= 80%)" if frT >= 0.8 else "outside (< 80%)",
    "supporting", f"TP1 r/mu(cal) {float(ft['r_hat']):.4g}{' AT GRID EDGE' if ft['edge'] == '1' else ''}; "
                  f"above/below {abT}/{beT}; misfit 20-100 kb {mT[0]:+.2f} half-bands")
if frT < 0.8:
    misfits.append("TP1 LD curve")
fsh = fit[("BAM", "primary", "in_sample_shape")]
add("BAM shape-only fit (23e, in-sample)", "genome", float(fsh["frac_in_band_own"]), "", "",
    "inside (>= 80%)" if fsh["adequate_own_80"] == "1" else "outside (< 80%)", "supporting",
    f"shape r/mu {float(fsh['r_hat']):.4g}{' AT GRID EDGE' if fsh['edge'] == '1' else ''}")
if fsh["adequate_own_80"] != "1":
    misfits.append("BAM shape-only fit")
# recent BAM Ne (last 10 generations), both mu arms
P = json.load(open(pj))
for mu in ("1e-08", "4e-09"):
    q = P.get(f"{model}|{mu}")
    if q is None:
        continue
    if "N_BP_2" in q:
        ne = [q["N_BP_2"]] * 10 if q.get("T_c_gen", 0) >= 10 else None
    elif model.rstrip("f") == "M3":
        gr = np.log(q["N_BP_end"] / q["N_BP_start"]) / q["T_gen"]
        ne = [q["N_BP_end"] * np.exp(-gr * t) for t in range(10)]
    else:
        ne = [q["N_BP_end"]] * 10
    hm = 10 / np.sum(1 / np.array(ne))
    add(f"recent BAM Ne (mu {mu})", "model", hm, 3100, 1500,
        ("above both temporal-Ne bounds" if hm >= 3100 else ("above g = 2 yr bound only" if hm >= 1500 else "below both")),
        "supporting", "harmonic mean over generations 0-9; temporal-Ne bounds >= ~3,100 (g = 1 yr), >= ~1,500 (g = 2 yr)")
# sensitivity rows and fixed r/mu points (scoring family)
fam_pr = A.family(fs, model, "ho", "primary", "BAM")
fr_pr_fam = A.bins_in_band(A.interp_at(fam_pr, r_hat), band["lo"], band["hi"])[0]
for t, row in (("callable", "callable"), ("gqdrop", "primary"), ("d019", "d019"), ("eps0.005", "primary"),
               ("eps0.01", "primary")):
    ftt = fit[("BAM", t, "calibration")]
    rt = float(ftt["r_hat"])
    bnd = obs[row][("BAM", "heldout")]
    frt = A.bins_in_band(A.interp_at(A.family(fs, model, "ho", t, "BAM"), rt), bnd["lo"], bnd["hi"])[0]
    lab = {"callable": "callable-only (sensitivity)", "gqdrop": "uncalibrated sensitivity (GQ-derived)",
           "d019": "strict-genotype mask (observed and simulated)", "eps0.005": "symmetric eps 0.005 (no empirical basis)",
           "eps0.01": "symmetric eps 0.01 (no empirical basis)"}[t]
    add(f"sensitivity: {lab}", "heldout", frt, "", "", f"bins-in-band shift {frt - fr_pr_fam:+.2f}", "never gating",
        f"r/mu(cal) {rt:.4g} [{float(ftt['lo95']):.4g}, {float(ftt['hi95']):.4g}]"
        f"{' AT GRID EDGE' if ftt['edge'] == '1' else ''}, shift vs primary {rt - r_hat:+.4g}; primary on the same "
        f"(scoring-family) basis {fr_pr_fam:.2f}", r_over_mu=rt, r_lo95=ftt["lo95"], r_hi95=ftt["hi95"], edge=ftt["edge"])
for rr in (0.25, 1.0, 2.5):
    c_ = A.family(fs, model, "ho", "primary", "BAM")[rr]
    f_, a_, b_ = A.bins_in_band(c_, band["lo"], band["hi"])
    add(f"P1 held-out at fixed point r/mu = {rr:g}", "heldout", f_, "", "", "80% rule " + ("met" if f_ >= 0.8 else
                                                                                               "not met"),
        "reported", f"above/below {a_}/{b_}; misfit 20-100 kb {misfit(c_, band)[0]:+.2f} half-bands; {COMP}",
        r_over_mu=rr)
passes = bool(p1 and admiss and (p3 or not p3g))
verdict = "passes" if passes else "fails"
if passes and misfits:
    verdict = "passes, with misfit in " + "; ".join(misfits)
fails = [r["criterion"] for r in rows if r["gating"] == "gating" and r["outcome"].startswith("fail")]
add("VERDICT", "-", "", "", "", verdict, "-", ("failing gating criteria: " + "; ".join(fails)) if fails else "")
cols = ["model", "analysis", "criterion", "set", "value", "lo", "hi", "outcome", "gating", "note", "r_over_mu", "r_lo95",
        "r_hi95", "edge"]
with open(outp + "_scores.tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", lineterminator="\n", restval="")
    w.writeheader(); w.writerows(rows)
np.savez(outp + "_curves.npz", sim_ho_BAM=sim_ho, sim_cal_BAM=pooled("BAM", cal_c), sim_ho_TP1=pooled("TP1", ho_c),
         sim_cal_TP1=pooled("TP1", cal_c),
         ho_callable_BAM=A.interp_at(A.family(fs, model, "ho", "callable", "BAM"),
                                     float(fit[("BAM", "callable", "calibration")]["r_hat"])))
print(model, verdict)
