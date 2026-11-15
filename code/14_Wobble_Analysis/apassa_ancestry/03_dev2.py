#!/usr/bin/env python3
"""03_dev2.py -- corroboration and shared-region tests with independent baselines.

Phased excess e = posterior fraction minus the B9 mean (leave-one-out for B9 fish); unphased excess x' = projection
score minus the P9 mean (P9 = the 9 Apassa 2009 TP fish that were local_ancestry primary-panel members; leave-one-out for P9
fish). Corroboration per fish (12 wobbles + 9 B9 controls): Spearman rho(e, x') over 2-Mb windows with a rotation
null, and the stratum-matched background test on x'. Shared regions: Delta S = mean pairwise Spearman of x' among
2019 wobbles minus that among B9 fish, null by random relabelling of the pooled fish (10,000; seed 32033).

Usage (repo root): python3 code/14_Wobble_Analysis/apassa_ancestry/03_dev2.py <outdir> <la_dir>
"""
import sys, csv, json, itertools
import numpy as np
from scipy.stats import rankdata, spearmanr

OUT, LA = sys.argv[1:3]
NPERM = 10_000
rng = np.random.default_rng(32033)
z = np.load(f"{OUT}/window_profiles.npz")
samp = {r["sample"]: r for r in csv.DictReader(open(f"{LA}/samples.tsv"), delimiter="\t")}
wins = z["wins"]; fP = list(z["fish_phased"]); fU = list(z["fish_unphased"])
B9 = list(z["B9"]); W19 = list(z["W19"]); W17 = list(z["W17"]); W09 = list(z["W09"])
P9 = [s for s in fU if samp[s]["group"] == "APA2009_TP" and s not in B9]
assert len(P9) == 9 and len(B9) == 9
Vw, Vc = z["V_phased_win"], z["V_phased_chrom"]
Qw, Qc = z["Q_unphased_win"], z["Q_unphased_chrom"]
ip = {f: k for k, f in enumerate(fP)}; iu = {f: k for k, f in enumerate(fU)}
ninf, ncall, nsite = z["ninf"], z["ncall"], z["nsite"]

def ex_phased(V, f):
    b = [x for x in B9 if x != f]
    return V[:, ip[f]] - V[:, [ip[x] for x in b]].mean(1)
def ex_unphased(Q, f):
    b = [x for x in P9 if x != f]
    return Q[:, iu[f]] - np.nanmean(Q[:, [iu[x] for x in b]], 1)

def zrank(v):
    r = rankdata(v); return (r - r.mean()) / r.std()
tert = lambda v: np.digitize(v, np.nanpercentile(v, [100 / 3, 200 / 3]))
RES = {"P9": P9}
rows = []
for f in W19 + W09 + B9:
    e = ex_phased(Vw, f); x = ex_unphased(Qw, f)
    ok = np.isfinite(e) & np.isfinite(x); e_, x_ = e[ok], x[ok]; U = int(ok.sum())
    rho = float(spearmanr(e_, x_)[0])
    re_, rx_ = zrank(e_), zrank(x_)
    rn = np.array([float((re_ * np.roll(rx_, k)).mean()) for k in rng.integers(1, U, NPERM)])
    p_rho = float((1 + (rn >= rho).sum()) / (NPERM + 1))
    k = iu[f]
    strata = (tert(ninf[ok, k]) * 9 + tert(ncall[ok, k] / np.maximum(ninf[ok, k], 1)) * 3 + tert(nsite[ok]))
    top = e_ >= np.percentile(e_, 90)
    def stat(tp):
        num = den = 0.0
        for st in np.unique(strata[tp]):
            inS = strata == st; nt = (tp & inS).sum(); bg = (~tp) & inS
            if bg.sum():
                num += nt * (x_[tp & inS].mean() - x_[bg].mean()); den += nt
        return num / den
    m0 = stat(top); mn = np.empty(NPERM)
    groups = [(np.flatnonzero(strata == st), int(top[strata == st].sum())) for st in np.unique(strata)]
    for p in range(NPERM):
        tp = np.zeros(U, bool)
        for ii, kk in groups:
            if kk: tp[rng.choice(ii, kk, replace=False)] = True
        mn[p] = stat(tp)
    p_bg = float((1 + (mn >= m0).sum()) / (NPERM + 1))
    ec = ex_phased(Vc, f); xc = ex_unphased(Qc, f)
    rows.append(dict(fish=f, group=samp[f]["group"], site=samp[f]["site"], role="control_B9" if f in B9 else "target",
                     n_windows=U, rho=rho, p_rho=p_rho, top_minus_matched_bg=float(m0), p_matched_bg=p_bg,
                     both_p_le_0_05=int(p_rho <= 0.05 and p_bg <= 0.05), rho_chrom=float(spearmanr(ec, xc)[0]),
                     excess_gw_phased=float(np.nanmean(e)), excess_gw_unphased_P9=float(np.nanmean(x))))
keys = list(rows[0])
with open(f"{OUT}/dev2_corroboration.tsv", "w") as fh:
    fh.write("\t".join(keys) + "\n")
    for r in rows:
        fh.write("\t".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
ctrl_max = max(r["rho"] for r in rows if r["role"] == "control_B9")
t19 = [r for r in rows if r["fish"] in W19]
RES["corroboration"] = dict(control_rho_range=[min(r["rho"] for r in rows if r["role"] == "control_B9"), ctrl_max],
                            control_both_sig=sum(r["both_p_le_0_05"] for r in rows if r["role"] == "control_B9"),
                            W19_both_sig=sum(r["both_p_le_0_05"] for r in t19),
                            W19_rho_above_control_max=sum(r["rho"] > ctrl_max for r in t19),
                            W17_rho_above_control_max=sum(r["rho"] > ctrl_max for r in t19 if r["fish"] in W17))
c = RES["corroboration"]
c["class"] = ("corroborated" if (c["W19_both_sig"] >= 8 and c["W19_rho_above_control_max"] >= 8) else
              "not_corroborated" if (c["W19_both_sig"] <= 2) else "partial")

# shared regions: Delta S with relabelling null
X = {f: ex_unphased(Qw, f) for f in W19 + B9 + W09}
okw = np.all(np.isfinite(np.stack([X[f] for f in W19 + B9], 1)), 1)
def meanS(fl, M):
    R = np.stack([zrank(M[f][okw]) for f in fl], 1); C = (R.T @ R) / R.shape[0]
    return float(C[np.triu_indices(len(fl), 1)].mean())
def dS_test(targets, controls, M, tag):
    pool = targets + controls; nt = len(targets)
    R = np.stack([zrank(M[f][okw]) for f in pool], 1); C = (R.T @ R) / R.shape[0]
    def ms(ix):
        sub = C[np.ix_(ix, ix)]; return sub[np.triu_indices(len(ix), 1)].mean()
    obs = ms(list(range(nt))) - ms(list(range(nt, len(pool))))
    null = np.empty(NPERM)
    for p in range(NPERM):
        perm = rng.permutation(len(pool)); null[p] = ms(list(perm[:nt])) - ms(list(perm[nt:]))
    RES[f"dS_{tag}"] = dict(S_targets=float(ms(list(range(nt)))), S_controls=float(ms(list(range(nt, len(pool))))),
                            S_cross=float(C[np.ix_(range(nt), range(nt, len(pool)))].mean()), dS=float(obs),
                            null_p95=float(np.percentile(null, 95)), p=float((1 + (null >= obs).sum()) / (NPERM + 1)),
                            n_windows=int(okw.sum()))
dS_test(W19, B9, X, "unphased_W19_vs_B9")
dS_test(W17, B9, X, "unphased_W17_vs_B9")
E = {f: ex_phased(Vw, f) for f in W19 + B9 + W09}
dS_test(W19, B9, E, "phased_W19_vs_B9_descriptive")
# windows exceeding the P9-based 95th percentile in >= 8/10 (>= 6/7)
P9x = np.stack([ex_unphased(Qw, f) for f in P9], 1); thr = np.nanpercentile(P9x, 95, axis=1)
lst = []
for tag, members, kmin in (("W19", W19, 8), ("W17", W17, 6), ("B9_control", B9, 7)):
    M = np.stack([X[f] for f in members], 1); n = (M > thr[:, None]).sum(1)
    hits = np.flatnonzero((n >= kmin) & okw)
    RES[f"n_windows_exceed_{tag}"] = int(hits.size)
    for h in hits:
        lst.append(dict(set=tag, window=str(wins[h]), n_exceed=int(n[h]), n_of=len(members), mean_x_unphased=float(np.nanmean(M[h])),
                        mean_e_phased=float(np.mean([E[f][h] for f in members]))))
with open(f"{OUT}/dev2_shared_windows.tsv", "w") as fh:
    fh.write("set\twindow\tn_exceed\tn_of\tmean_x_unphased\tmean_e_phased\n")
    for r in lst:
        fh.write(f"{r['set']}\t{r['window']}\t{r['n_exceed']}\t{r['n_of']}\t{r['mean_x_unphased']:.4f}\t{r['mean_e_phased']:.4f}\n")
RES["shared_class"] = "recurrent_shared_regions" if RES["dS_unphased_W19_vs_B9"]["p"] <= 0.05 else "no_excess_sharing"
json.dump(RES, open(f"{OUT}/dev2_results.json", "w"), indent=1)
print(json.dumps(RES, indent=1))
