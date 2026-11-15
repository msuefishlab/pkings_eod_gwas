#!/usr/bin/env python3
"""03_markers.py -- ancestry-informative markers per reference panel.

From the candidate sites of 02 (phased haplotypes of all 114 APA/BAM/BEN samples):
  1. AIM: |p_BAM - p_ASC| >= 0.5, where ASC is an ascertainment TP sample disjoint from the panel's TP emission
     sample (picking AIMs on the emission panel itself gives a winner's curse).
  2. Thin: one AIM per non-overlapping 10-kb bin, the one with the largest |delta| (lowest position on ties).
  3. Prune: walking along each chromosome, drop an AIM whose haplotype r^2 with the previous kept AIM exceeds
     0.5 in either panel.
Each marker is flagged gw (outside the six extended peaks +/- 250 kb and the chr4 inversion, used for fitting
and genome-wide summaries) and with the peak id whose extended interval contains it (0 if none).

Emission frequencies (la_lib.load_markers) are empirical-Bayes: per panel and side, a beta prior fitted by
maximum beta-binomial likelihood to the counts at the kept genome-wide markers, oriented to the Bambomo-high
allele; the prior parameters are stored here.

Panels (BAM emission, TP emission, TP ascertainment): primary (PB, PA, BN), full (FB, FA, BN), ben (PB, BN, FA);
membership from samples.tsv (01_panels.py).
Writes <out_dir>/markers_<panel>.npz: chrom, pos, H (markers x haplotypes, 0/1), hap_sample, cB, nB, cT, nT,
orient (1 if allele 1 is the Bambomo-high allele at ascertainment), prior_B, prior_T (alpha, beta), gw, peak;
and markers_summary.tsv.

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/03_markers.py <out_dir>
"""
import gzip, sys
import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln

PANELS = {"primary": ("PB", "PA", "BN"), "full": ("FB", "FA", "BN"), "ben": ("PB", "BN", "FA")}


def beta_prior(c, n):
    """Maximum beta-binomial likelihood (alpha, beta) for counts c of n."""
    nll = lambda x: -(betaln(c + np.exp(x[0]), n - c + np.exp(x[1])) - betaln(np.exp(x[0]), np.exp(x[1]))).sum()
    r = minimize(nll, [0.0, 0.0], method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-6, maxiter=2000))
    return np.exp(r.x)
INV = ("chr4", 40_880_000, 42_770_000)
BIN, DMIN, R2MAX = 10_000, 0.5, 0.5


def regions():
    excl, peaks = [INV], []
    for line in open("input_data/10_Introgression/peaks.bed"):
        f = line.split()
        ch, a, b, pid = f"chr{f[0]}", int(f[1]), int(f[2]), int(f[3])
        excl.append((ch, a - 250_000, b + 250_000))
        peaks.append((ch, a, b, pid))
    return excl, peaks


def r2(x, y):
    x = x.astype(float); y = y.astype(float)
    vx, vy = x.var(), y.var()
    if vx == 0 or vy == 0:
        return 0.0
    c = ((x - x.mean()) * (y - y.mean())).mean()
    return c * c / (vx * vy)


def main():
    out = sys.argv[1]
    samples = open(f"{out}/cand/samples.txt").read().split()
    roles = {}
    for line in open(f"{out}/samples.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "sample":
            roles[f[0]] = set(f[5].split(","))
    hap_sample = np.repeat(np.array(samples), 2)
    hcols = {g: np.flatnonzero([g in roles.get(s, ()) for s in hap_sample]) for g in ("PB", "PA", "FB", "FA", "BN")}
    excl, peaks = regions()
    data = {}
    for c in range(1, 26):
        ch = f"chr{c}"
        pos, rows = [], []
        with gzip.open(f"{out}/cand/{ch}.tsv.gz", "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                gt = f[2:]
                assert all(len(g) == 3 and g[1] == "|" for g in gt), (ch, f[1])
                pos.append(int(f[1]))
                rows.append([int(a) for g in gt for a in (g[0], g[2])])
        data[ch] = (np.array(pos), np.array(rows, dtype=np.int8))
        print(ch, len(pos), flush=True)
    summ = ["panel\tchrom\tcandidates\taims\tthinned\tpruned\tgw\tpeak\tmean_adj_r2_B\tmean_adj_r2_T"]
    for panel, (gb, gt_, ga) in PANELS.items():
        cb, ct, ca = hcols[gb], hcols[gt_], hcols[ga]
        keep_all = {k: [] for k in ("chrom", "pos", "H", "cB", "nB", "cT", "nT", "orient", "gw", "peak")}
        for c in range(1, 26):
            ch = f"chr{c}"
            pos, H = data[ch]
            pB = H[:, cb].mean(1); pA = H[:, ca].mean(1)
            d = np.abs(pB - pA)
            aim = np.flatnonzero(d >= DMIN)
            # thin: best per 10-kb bin (lowest position on ties: stable sort on -d within bin)
            b = (pos[aim] - 1) // BIN
            order = np.lexsort((pos[aim], -d[aim], b))
            first = np.ones(order.size, bool); first[1:] = b[order][1:] != b[order][:-1]
            thin = np.sort(aim[order[first]])
            kept, rB, rT = [], [], []
            for j in thin:
                if kept:
                    k = kept[-1]
                    a, bb = r2(H[j, cb], H[k, cb]), r2(H[j, ct], H[k, ct])
                    if a > R2MAX or bb > R2MAX:
                        continue
                    rB.append(a); rT.append(bb)
                kept.append(j)
            kept = np.array(kept, int)
            p = pos[kept]
            gw = np.ones(p.size, bool)
            for x, a, bb in excl:
                if x == ch:
                    gw &= ~((p >= a) & (p <= bb))
            pk = np.zeros(p.size, int)
            for x, a, bb, pid in peaks:
                if x == ch:
                    pk[(p >= a) & (p <= bb)] = pid
            keep_all["chrom"].append(np.full(p.size, c, np.int8)); keep_all["pos"].append(p)
            keep_all["H"].append(H[kept])
            keep_all["cB"].append(H[kept][:, cb].sum(1)); keep_all["nB"].append(np.full(p.size, cb.size))
            keep_all["cT"].append(H[kept][:, ct].sum(1)); keep_all["nT"].append(np.full(p.size, ct.size))
            keep_all["orient"].append((pB[kept] >= pA[kept]).astype(np.int8))
            keep_all["gw"].append(gw); keep_all["peak"].append(pk)
            summ.append(f"{panel}\t{ch}\t{pos.size}\t{aim.size}\t{thin.size}\t{kept.size}\t{gw.sum()}\t{(pk > 0).sum()}\t"
                        f"{np.mean(rB) if rB else float('nan'):.4f}\t{np.mean(rT) if rT else float('nan'):.4f}")
        arr = {k: np.concatenate(v) for k, v in keep_all.items()}
        o, g = arr["orient"].astype(bool), arr["gw"].astype(bool)
        ob = np.where(o, arr["cB"], arr["nB"] - arr["cB"]); ot = np.where(o, arr["cT"], arr["nT"] - arr["cT"])
        arr["prior_B"] = beta_prior(ob[g], arr["nB"][g]); arr["prior_T"] = beta_prior(ot[g], arr["nT"][g])
        print(panel, "beta priors (Bambomo-high allele): BAM", arr["prior_B"].round(3), "TP", arr["prior_T"].round(3))
        np.savez_compressed(f"{out}/markers_{panel}.npz", hap_sample=hap_sample, **arr)
        tot = [l.split("\t") for l in summ[1:] if l.startswith(panel + "\t")]
        print(panel, "aims", sum(int(t[3]) for t in tot), "thinned", sum(int(t[4]) for t in tot),
              "pruned", sum(int(t[5]) for t in tot), "gw", sum(int(t[6]) for t in tot),
              "peak", sum(int(t[7]) for t in tot), flush=True)
    open(f"{out}/markers_summary.tsv", "w").write("\n".join(summ) + "\n")


if __name__ == "__main__":
    main()
