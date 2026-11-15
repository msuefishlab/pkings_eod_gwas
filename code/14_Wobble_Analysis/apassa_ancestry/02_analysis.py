#!/usr/bin/env python3
"""02_analysis.py -- Parts A and B of the Apassa ancestry-distribution analysis. numpy/scipy only.

Part A (phased, local_ancestry primary-panel outputs): per-fish ancestry fraction, tract statistics, excess over the
2009 Apassa TP baseline (B9; B4 sensitivity) per chromosome and 2-Mb window, concentration shares, the dispersion
statistic D against local_ancestry-simulator fish at each fish's matched cell (01_sim_null.py), and shared-window
tests among the 2019 wobbles (circular-rotation null).
Part B (phase-free): unphased genotypes from simulation/26a_extract_temporal.sh (DP >= 8, GQ >= 20, no imputation);
per-fish projection scores per chromosome and window against Bambomo 2009 BP and Apassa 2009 TP (leave-one-out for
TP targets); gates; shared-window tests; corroboration of the phased excess (rotation and stratum-matched background).

Usage (repo root): python3 code/14_Wobble_Analysis/apassa_ancestry/02_analysis.py <la_dir> <temporal_dir> <qc_tsv> <outdir>
"""
import os, sys, csv, gzip, json, itertools
import numpy as np
from scipy.stats import rankdata, spearmanr, pearsonr

sys.path.insert(0, "code/11_Drift_vs_Selection/local_ancestry")
import la_lib as L

LA, D19, QC, OUT = sys.argv[1:5]
NPERM, SEED, WIN, CAP = 10_000, 32032, 2_000_000, 100_000
rng = np.random.default_rng(SEED)
RES = {}

def tsv(path, rows, keys=None):
    keys = keys or list(rows[0])
    with open(path, "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[k]:.6g}" if isinstance(r[k], (float, np.floating)) else str(r[k]) for k in keys) + "\n")

samp = {r["sample"]: r for r in csv.DictReader(open(f"{LA}/samples.tsv"), delimiter="\t")}
pf = {r["sample"]: r for r in csv.DictReader(open(f"{LA}/results/per_fish.tsv"), delimiter="\t")}
grp = lambda s: samp[s]["group"]
site = lambda s: samp[s]["site"]

# ------------------------------------------------------------------ markers, callable weights, units
z = np.load(f"{LA}/hap_primary/post_gw.npz")
post, chrom, pos, hs = z["post"], z["chrom"].astype(int), z["pos"], z["hap_sample"]
fish = list(dict.fromkeys(hs.tolist()))
dose = np.stack([post[:, hs == f].sum(1) for f in fish], 1)          # (M, F) 0..2
seg = L.segments(chrom, pos, L.excluded_regions())
Mn = pos.size
w = np.zeros(Mn)
for i in range(Mn):
    left = min((pos[i] - pos[i - 1]) / 2, CAP) if i and seg[i - 1] == seg[i] else 0.0
    right = min((pos[i + 1] - pos[i]) / 2, CAP) if i + 1 < Mn and seg[i + 1] == seg[i] else 0.0
    w[i] = left + right
win_id = np.array([f"chr{c}:{p // WIN}" for c, p in zip(chrom, pos)])
CH = np.arange(1, 26)
wins = list(dict.fromkeys(win_id.tolist()))
win_idx = {u: k for k, u in enumerate(wins)}
wi = np.array([win_idx[u] for u in win_id])
def unit_vals(D, uidx, nu):
    num = np.zeros((nu, D.shape[1])); den = np.zeros(nu)
    np.add.at(num, uidx, w[:, None] * D / 2); np.add.at(den, uidx, w)
    return num / den[:, None], den / den.sum()
ci = chrom - 1
Vc, wc = unit_vals(dose, ci, 25)
Vw, ww = unit_vals(dose, wi, len(wins))
RES["callable_Mb"] = float(w.sum() / 1e6); RES["analysed_Mb"] = L.analysed_length(pos, seg) / 1e6
RES["n_markers"] = int(Mn); RES["n_windows"] = len(wins)

APA = [f for f in fish if grp(f).startswith("APA")]
W19 = [f for f in fish if grp(f) == "APA2019_Wobble"]; W17 = [f for f in W19 if site(f) == "2019-017"]
W09 = [f for f in fish if grp(f) == "APA2009_Wobble"]
B9 = [f for f in fish if grp(f) == "APA2009_TP"]
B4 = [f for f in B9 if "holdout_APA" in samp[f]["roles"]]
assert len(W19) == 10 and len(W09) == 2 and len(B9) == 9 and len(B4) == 4, (len(W19), len(W09), len(B9), len(B4))
fi = {f: k for k, f in enumerate(fish)}

def excess(V, base, f):
    b = [x for x in base if x != f]
    return V[:, fi[f]] - V[:, [fi[x] for x in b]].mean(1)

# ------------------------------------------------------------------ Part A per-fish table
tracts = {}
for r in csv.DictReader(open(f"{LA}/hap_primary/tracts.tsv"), delimiter="\t"):
    tracts.setdefault(r["sample"], []).append((int(str(r["chrom"]).replace("chr", "")), float(r["start"]), float(r["end"]), float(r["length"])))
alen = L.analysed_length(pos, seg)
B9loo_c = np.stack([excess(Vc, B9, b) for b in B9], 1)              # (25, 9)
thr_c = np.percentile(B9loo_c, 95, axis=1)
B9loo_w = np.stack([excess(Vw, B9, b) for b in B9], 1)
thr_w = np.percentile(B9loo_w, 95, axis=1)

def topk_share(vals, wts, k):
    pos_ex = np.maximum(vals, 0) * wts
    tot = pos_ex.sum()
    return float(np.sort(pos_ex)[::-1][:k].sum() / tot) if tot > 0 else np.nan

def disp(V, wts, k):
    v = V[:, k]; vb = (wts * v).sum()
    return float((wts * (v - vb) ** 2).sum())

def tract_stats(tl):
    tl = np.sort(np.asarray(tl))[::-1]
    tot = tl.sum()
    return dict(n_tracts=int(tl.size), tract_median_kb=float(np.median(tl) / 1e3) if tl.size else np.nan,
                tract_mean_kb=float(tl.mean() / 1e3) if tl.size else np.nan, tract_max_Mb=float(tl[0] / 1e6) if tl.size else 0.0,
                top5_tract_share=float(tl[:5].sum() / tot) if tot else np.nan, top10_tract_share=float(tl[:10].sum() / tot) if tot else np.nan,
                n_tracts_ge2Mb=int((tl >= 2e6).sum()), tract_cov=float(tot / (2 * alen)),
                tract_cov_ge2Mb=float(tl[tl >= 2e6].sum() / (2 * alen)))

rowsA = []
for f in APA:
    k = fi[f]
    ec = excess(Vc, B9, f); ew = excess(Vw, B9, f)
    ec4 = excess(Vc, B4, f)
    tl = [t[3] for t in tracts.get(f, [])]
    r = dict(fish=f, group=grp(f), site=site(f), q_m=float(pf[f]["q_m"]), match_T=pf[f]["match_T"], match_m=pf[f]["match_m"],
             post_frac_marker=float(dose[:, k].mean() / 2), post_frac_callable=float((w * dose[:, k]).sum() / (2 * w.sum())),
             **tract_stats(tl),
             excess_gw_B9=float((wc * ec).sum()), excess_gw_B4=float((wc * ec4).sum()),
             n_chrom_detectable=int((ec > thr_c).sum()), n_chrom_with_ge2Mb_tract=len({t[0] for t in tracts.get(f, []) if t[3] >= 2e6}),
             E_pos=float((np.maximum(ec, 0) * wc).sum()),
             top1_chrom_share=topk_share(ec, wc, 1), top3_chrom_share=topk_share(ec, wc, 3), top5_chrom_share=topk_share(ec, wc, 5),
             top5_win_share=topk_share(ew, ww, 5), top10_win_share=topk_share(ew, ww, 10),
             D_chrom=disp(Vc, wc, k), D_win=disp(Vw, ww, k))
    rowsA.append(r)

# ------------------------------------------------------------------ simulation null (01_sim_null.py)
zs = np.load(f"{OUT}/sim_null.npz")
assert np.array_equal(zs["pos"], pos) and np.array_equal(zs["chrom"].astype(int), chrom)
sd = zs["dose"].astype(float); cT, cM = zs["cell_T"], zs["cell_m"]
SVc, _ = unit_vals(sd, ci, 25); SVw, _ = unit_vals(sd, wi, len(wins))
base0 = np.flatnonzero(cM == 0)
S_ec = SVc - SVc[:, base0].mean(1, keepdims=True); S_ew = SVw - SVw[:, base0].mean(1, keepdims=True)
sim_rows = []
for j in range(sd.shape[1]):
    tl = []
    for cal in (zs["call0"][:, j], zs["call1"][:, j]):
        for s_, a, b_ in L.call_tracts(cal.astype(float), pos, seg):
            tl.append(b_ - a)
    sim_rows.append(dict(T=cT[j], m=cM[j], D_chrom=disp(SVc, wc, j), D_win=disp(SVw, ww, j),
                         top5_chrom_share=topk_share(S_ec[:, j], wc, 5), top5_win_share=topk_share(S_ew[:, j], ww, 5),
                         excess_gw=float((wc * S_ec[:, j]).sum()), **tract_stats(tl)))
tsv(f"{OUT}/sim_null_fish.tsv", sim_rows)
for r in rowsA:
    T_, m_ = float(r["match_T"]), float(r["match_m"])
    if m_ == 0: T_ = 0.0
    idx = [j for j, s in enumerate(sim_rows) if s["T"] == T_ and abs(s["m"] - m_) < 1e-9]
    r["n_sim_matched"] = len(idx)
    for key in ("D_chrom", "D_win", "top5_chrom_share", "top5_win_share", "top5_tract_share", "tract_max_Mb", "n_tracts_ge2Mb"):
        sv = np.array([sim_rows[j][key] for j in idx], float); sv = sv[np.isfinite(sv)]
        r[f"p_{key}"] = float((sv >= r[key]).mean()) if sv.size and np.isfinite(r[key]) else np.nan
        r[f"sim_median_{key}"] = float(np.median(sv)) if sv.size else np.nan
    pc, pw = r["p_D_chrom"], r["p_D_win"]
    r["dist_class"] = ("concentrated" if (pc <= 0.01 or pw <= 0.01) else
                       "consistent_with_random_tracts" if (pc > 0.05 and pw > 0.05) else "intermediate")
tsv(f"{OUT}/partA_per_fish.tsv", rowsA)
for nm, grpf in (("W19", W19), ("W17", W17)):
    cls = [r["dist_class"] for r in rowsA if r["fish"] in grpf]
    RES[f"dist_class_counts_{nm}"] = {c: cls.count(c) for c in set(cls)}
# per-chromosome excess table
rows = []
for f in APA:
    ec = excess(Vc, B9, f)
    for c in range(25):
        rows.append(dict(fish=f, group=grp(f), site=site(f), chrom=c + 1, v=float(Vc[c, fi[f]]), excess=float(ec[c]),
                         B9_p95=float(thr_c[c]), detectable=int(ec[c] > thr_c[c]), callable_Mb=float(wc[c] * w.sum() / 1e6)))
tsv(f"{OUT}/partA_chrom_excess.tsv", rows)

# ------------------------------------------------------------------ shared-window tests (rotation null)
def shared_tests(E, thr, members, tag, kmin):
    """E: (U, n) excess profiles; thr: (U,) per-window B9 95th pct; rotation null on each fish independently."""
    U, n = E.shape
    R = np.stack([rankdata(E[:, j]) for j in range(n)], 1)
    R = (R - R.mean(0)) / R.std(0)
    X = E > thr[:, None]
    def S_of(Rm):
        C = (Rm.T @ Rm) / U; iu = np.triu_indices(n, 1); return float(C[iu].mean())
    S0 = S_of(R); K0 = int((X.sum(1) >= kmin).sum())
    Sn = np.empty(NPERM); Kn = np.empty(NPERM, int)
    for p in range(NPERM):
        sh = rng.integers(0, U, n)
        Rr = np.stack([np.roll(R[:, j], sh[j]) for j in range(n)], 1)
        Xr = np.stack([np.roll(X[:, j], sh[j]) for j in range(n)], 1)
        Sn[p] = S_of(Rr); Kn[p] = int((Xr.sum(1) >= kmin).sum())
    out = dict(S=S0, S_null_median=float(np.median(Sn)), S_null_p95=float(np.percentile(Sn, 95)),
               p_S=float((1 + (Sn >= S0).sum()) / (NPERM + 1)), K=K0, kmin=kmin, K_null_median=float(np.median(Kn)),
               K_null_p95=float(np.percentile(Kn, 95)), p_K=float((1 + (Kn >= K0).sum()) / (NPERM + 1)),
               n_windows=U, n_fish=n)
    hits = np.flatnonzero(X.sum(1) >= kmin)
    RES[f"shared_{tag}"] = out
    return out, hits, X.sum(1)

EW = {f: excess(Vw, B9, f) for f in APA}
shared_rows = []
for tag, members, kmin in (("phased_W19", W19, 8), ("phased_W17", W17, 6)):
    E = np.stack([EW[f] for f in members], 1)
    out, hits, nexc = shared_tests(E, thr_w, members, tag, kmin)
    for h in hits:
        shared_rows.append(dict(test=tag, window=wins[h], n_exceed=int(nexc[h]), mean_excess=float(E[h].mean())))
# also the windows with the highest carrier counts (descriptive, top 15 per test)
for tag, members in (("phased_W19_top", W19), ("phased_W17_top", W17)):
    E = np.stack([EW[f] for f in members], 1); nexc = (E > thr_w[:, None]).sum(1)
    for h in np.argsort(-nexc - E.mean(1) * 1e-3)[:15]:
        shared_rows.append(dict(test=tag, window=wins[h], n_exceed=int(nexc[h]), mean_excess=float(E[h].mean())))

# ------------------------------------------------------------------ Part B: unphased projection scores
order = [l.strip() for l in open(f"{D19}/samples.txt") if l.strip()]
col = {s: k for k, s in enumerate(order)}
BBP = [s for s in order if samp.get(s, {}).get("group") == "BAM2009_BP"]
ATP = [s for s in order if samp.get(s, {}).get("group") == "APA2009_TP"]
TARG = [s for s in order if samp.get(s, {}).get("group", "").startswith("APA")]
assert len(BBP) == 41 and len(ATP) == 18 and len(TARG) == 30, (len(BBP), len(ATP), len(TARG))
iB = np.array([col[s] for s in BBP]); iA = np.array([col[s] for s in ATP])
numC = {s: np.zeros(25) for s in TARG}; denC = {s: np.zeros(25) for s in TARG}
numW = {s: np.zeros(len(wins)) for s in TARG}; denW = {s: np.zeros(len(wins)) for s in TARG}
ninfW = {s: np.zeros(len(wins)) for s in TARG}; ncallW = {s: np.zeros(len(wins)) for s in TARG}
nsiteW = np.zeros(len(wins)); ninf_gw = {}
hetinf = {s: [0, 0] for s in TARG}
GT = {b"0/0": 0, b"0/1": 1, b"1/0": 1, b"1/1": 2, b"0|0": 0, b"0|1": 1, b"1|0": 1, b"1|1": 2}
for c in range(1, 26):
    P, G = [], []
    with gzip.open(f"{D19}/geno/chr{c}.tsv.gz", "rb") as fh:
        for line in fh:
            x = line.rstrip(b"\n").split(b"\t")
            P.append(int(x[1])); G.append([GT.get(g, -1) for g in x[2:]])
    P = np.array(P); G = np.array(G, np.int8); assert G.shape[1] == len(order)
    Gf = np.where(G < 0, np.nan, G.astype(float))
    wid = np.array([win_idx.get(f"chr{c}:{p // WIN}", -1) for p in P])
    okw = wid >= 0
    np.add.at(nsiteW, wid[okw], 1)
    def fr(ix):
        sub = Gf[:, ix]; n = np.isfinite(sub).sum(1)
        return np.nansum(sub, 1) / (2 * np.maximum(n, 1)), n
    pB, nB = fr(iB); pA_all, nA_all = fr(iA)
    okB = nB >= 0.5 * len(iB)
    for s in TARG:
        if s in ATP:
            ix = np.array([col[a] for a in ATP if a != s]); pA, nA = fr(ix); okA = nA >= 0.5 * ix.size
        else:
            pA, nA, okA = pA_all, nA_all, nA_all >= 0.5 * len(iA)
        d = pB - pA
        inf = okB & okA & (np.abs(d) >= 0.3) & okw
        g = Gf[:, col[s]]; called = inf & np.isfinite(g)
        num = np.where(called, (g / 2 - pA) * d, 0.0); den = np.where(called, d * d, 0.0)
        numC[s][c - 1] += num.sum(); denC[s][c - 1] += den.sum()
        np.add.at(numW[s], wid[called], num[called]); np.add.at(denW[s], wid[called], den[called])
        np.add.at(ninfW[s], wid[inf], 1); np.add.at(ncallW[s], wid[called], 1)
        hetinf[s][0] += int((called & (g == 1)).sum()); hetinf[s][1] += int(called.sum())
    print(f"chr{c}: {P.size} sites", flush=True)
qC = {s: numC[s] / denC[s] for s in TARG}
qW = {s: np.where(denW[s] > 0, numW[s] / np.where(denW[s] > 0, denW[s], 1), np.nan) for s in TARG}
qG = {s: float(numC[s].sum() / denC[s].sum()) for s in TARG}
allden = np.concatenate([denW[s] for s in TARG]); thr_info = float(np.percentile(allden[allden > 0], 5))
for s in TARG:
    qW[s] = np.where(denW[s] >= thr_info, qW[s], np.nan)
RES["partB_info_threshold"] = thr_info
# gates
common = [s for s in TARG if s in pf]
r_gate = pearsonr([qG[s] for s in common], [float(pf[s]["q_m"]) for s in common])[0]
b9g = [qG[s] for s in B9]
RES["partB_gate"] = dict(r_unphased_vs_q_m=float(r_gate), n=len(common), B9_q_range=[float(min(b9g)), float(max(b9g))],
                         pass_r=bool(r_gate >= 0.8), pass_B9=bool(max(abs(x) for x in b9g) <= 0.05))
RES["partB_gate"]["pass"] = RES["partB_gate"]["pass_r"] and RES["partB_gate"]["pass_B9"]
# depth check
dp = {r["sample"]: float(r["mean_dp"]) for r in csv.DictReader(open(QC), delimiter="\t") if r["label"] == "genome"}
het = {s: hetinf[s][0] / max(hetinf[s][1], 1) for s in TARG}
for nm, ss in (("APA2009_TP", ATP), ("APA2019_W", W19)):
    RES[f"partB_het_vs_depth_{nm}"] = float(spearmanr([het[s] for s in ss], [dp[s] for s in ss])[0])
# excess vs B9 (leave-one-out), per window and chromosome
def xex(Q, s):
    b = [x for x in B9 if x != s]
    return Q[s] - np.nanmean(np.stack([Q[x] for x in b], 1), 1)
rowsB = []
for s in TARG:
    xc = xex(qC, s)
    rowsB.append(dict(fish=s, group=grp(s), site=site(s), depth=dp.get(s, np.nan), q_unphased_gw=qG[s],
                      q_m_D020=float(pf[s]["q_m"]) if s in pf else np.nan, het_at_informative=het[s],
                      informative_called=int(sum(ncallW[s])), excess_gw_unphased=float(qG[s] - np.mean([qG[b] for b in B9 if b != s])),
                      top5_chrom_share_unphased=topk_share(xc, wc, 5)))
tsv(f"{OUT}/partB_per_fish.tsv", rowsB)
rows = []
for s in TARG:
    xc = xex(qC, s)
    for c in range(25):
        rows.append(dict(fish=s, group=grp(s), chrom=c + 1, q=float(qC[s][c]), excess=float(xc[c]), info=float(denC[s][c])))
tsv(f"{OUT}/partB_chrom_excess.tsv", rows)
XW = {s: xex(qW, s) for s in TARG}
B9x = np.stack([XW[b] for b in B9], 1); thr_xw = np.nanpercentile(B9x, 95, axis=1)
okwin = np.all(np.isfinite(np.stack([XW[s] for s in W19 + B9], 1)), 1) & np.isfinite(thr_xw)
RES["partB_windows_used"] = int(okwin.sum())
for tag, members, kmin in (("unphased_W19", W19, 8), ("unphased_W17", W17, 6)):
    E = np.stack([XW[s][okwin] for s in members], 1)
    out, hits, nexc = shared_tests(E, thr_xw[okwin], members, tag, kmin)
    wl = np.array(wins)[okwin]
    for h in hits:
        shared_rows.append(dict(test=tag, window=wl[h], n_exceed=int(nexc[h]), mean_excess=float(E[h].mean())))
tsv(f"{OUT}/shared_windows.tsv", shared_rows, ["test", "window", "n_exceed", "mean_excess"])

# ------------------------------------------------------------------ corroboration (phased e vs unphased x)
tert = lambda v: np.digitize(v, np.nanpercentile(v, [100 / 3, 200 / 3]))
corr_rows = []
for s in W19 + W09:
    e = EW[s]; x = XW[s]
    ok = np.isfinite(e) & np.isfinite(x)
    e_, x_ = e[ok], x[ok]; U = ok.sum()
    rho0 = spearmanr(e_, x_)[0]
    re_, rx_ = rankdata(e_), rankdata(x_); re_ = (re_ - re_.mean()) / re_.std(); rx_ = (rx_ - rx_.mean()) / rx_.std()
    sh = rng.integers(1, U, NPERM)
    rn = np.array([float((re_ * np.roll(rx_, k)).mean()) for k in sh])
    p_rho = float((1 + (rn >= rho0).sum()) / (NPERM + 1))
    top = e_ >= np.percentile(e_, 90)
    strata = (tert(ninfW[s][ok]) * 9 + tert(ncallW[s][ok] / np.maximum(ninfW[s][ok], 1)) * 3 + tert(nsiteW[ok]))
    def stat(tp):
        num = den = 0.0
        for st in np.unique(strata[tp]):
            inS = strata == st; nt = (tp & inS).sum(); bg = (~tp) & inS
            if bg.sum() == 0: continue
            num += nt * (x_[tp & inS].mean() - x_[bg].mean()); den += nt
        return num / den if den else np.nan
    m0 = stat(top)
    mn = np.empty(NPERM)
    for p in range(NPERM):
        tp = np.zeros(U, bool)
        for st in np.unique(strata):
            ii = np.flatnonzero(strata == st); k = int(top[ii].sum())
            if k: tp[rng.choice(ii, k, replace=False)] = True
        mn[p] = stat(tp)
    p_bg = float((1 + (mn >= m0).sum()) / (NPERM + 1))
    corr_rows.append(dict(fish=s, group=grp(s), site=site(s), n_windows=int(U), rho=float(rho0), p_rho=p_rho,
                          top_minus_matched_bg=float(m0), p_matched_bg=p_bg,
                          corroborated=int(p_rho <= 0.05 and p_bg <= 0.05)))
tsv(f"{OUT}/partB_corroboration.tsv", corr_rows)
n_cor = sum(r["corroborated"] for r in corr_rows if r["fish"] in W19)
n_none = sum(int(r["p_rho"] > 0.05 and r["p_matched_bg"] > 0.05) for r in corr_rows if r["fish"] in W19)
RES["corroboration_W19"] = dict(n_corroborated=n_cor, n_neither=n_none,
                                 cls="corroborated" if n_cor >= 8 else "not_corroborated" if n_none >= 8 else "partial")
n_cor7 = sum(r["corroborated"] for r in corr_rows if r["fish"] in W17)
RES["corroboration_W17"] = dict(n_corroborated=n_cor7, n=len(W17))
# chromosome-level agreement phased vs unphased, per fish (descriptive)
RES["chrom_level_rho"] = {s: float(spearmanr(excess(Vc, B9, s), xex(qC, s))[0]) for s in W19 + W09}
np.savez_compressed(f"{OUT}/window_profiles.npz", wins=np.array(wins), w_win=ww, fish_phased=np.array(APA),
                    E_phased=np.stack([EW[f] for f in APA], 1), thr_phased=thr_w,
                    fish_unphased=np.array(TARG), X_unphased=np.stack([XW[s] for s in TARG], 1), thr_unphased=thr_xw,
                    ninf=np.stack([ninfW[s] for s in TARG], 1), nsite=nsiteW,
                    V_phased_win=Vw[:, [fi[f] for f in APA]], V_phased_chrom=Vc[:, [fi[f] for f in APA]], w_chrom=wc,
                    Q_unphased_win=np.stack([qW[s] for s in TARG], 1), Q_unphased_chrom=np.stack([qC[s] for s in TARG], 1),
                    ncall=np.stack([ncallW[s] for s in TARG], 1), B9=np.array(B9), W19=np.array(W19), W17=np.array(W17),
                    W09=np.array(W09))
json.dump(RES, open(f"{OUT}/results.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print(json.dumps(RES, indent=1, default=str))
