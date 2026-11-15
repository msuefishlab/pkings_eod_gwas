#!/usr/bin/env python3
"""23g_unlinked_floor.py -- is the r^2 floor raised by structure or relatedness? (LD-decay diagnostic)

The LD fit of 23e failed: observed r^2 is flatter with distance than any simulated curve. Raised r^2 between
unlinked SNPs would point to structure or relatedness. Here: one 23a segment per chromosome (seed 1), its first
200 kb, up to 150 SNPs passing the ld23_lib filters; r^2 for every pair of SNPs on different chromosomes, jointly
called genotypes. Null: the same with individuals permuted within the second block (keeps n, missingness and
within-block LD; breaks any cross-block association).

Usage (repo root): 23g_unlinked_floor.py <vcf> <samples.txt> <23a_observed_ld_<pop>.tsv>  (prints one summary)
"""
import csv, importlib.util, os, sys
import numpy as np
sys.path.insert(0, "code/11_Drift_vs_Selection/simulation")
import ld23_lib as L
spec = importlib.util.spec_from_file_location("a", "code/11_Drift_vs_Selection/simulation/23a_observed_ld_decay.py")
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)

vcf, samples, obs_tsv = sys.argv[1:4]
blocks = []
for r in csv.DictReader(open(obs_tsv), delimiter="\t"):
    k = (r["chrom"], int(r["block"]))
    if k not in blocks:
        blocks.append(k)
rng = np.random.default_rng(1)
pick = {}
for ch, b in blocks:
    pick.setdefault(ch, []).append(b)
sel = [(ch, rng.choice(bs)) for ch, bs in sorted(pick.items())]       # one block per chromosome
data = []
for ch, b in sel:
    s = b * 1_000_000
    pos, G = A.genotypes(vcf, samples, ch, s, s + 200_000)
    M = (G >= 0).astype(np.int8); X = np.where(G >= 0, G, 0).astype(np.int8)
    keep = L.filter_snps(X, M)
    idx = np.flatnonzero(keep)
    idx = rng.choice(idx, min(150, idx.size), replace=False)
    data.append((ch, X[idx].astype(float), M[idx].astype(float)))

def cross(XA, MA, XB, MB):
    n = MA @ MB.T; sx = XA @ MB.T; sy = MA @ XB.T
    sxx = (XA ** 2) @ MB.T; syy = MA @ (XB ** 2).T; sxy = XA @ XB.T
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = sxy / n - sx * sy / n ** 2; vx = sxx / n - (sx / n) ** 2; vy = syy / n - (sy / n) ** 2
        r2 = cov ** 2 / (vx * vy)
    ok = (n >= L.MIN_JOINT) & np.isfinite(r2) & (vx > 0) & (vy > 0)
    return r2[ok]

obs, nul = [], []
for i in range(len(data)):
    for j in range(i + 1, len(data)):
        _, XA, MA = data[i]; _, XB, MB = data[j]
        obs.append(cross(XA, MA, XB, MB))
        p = rng.permutation(XB.shape[1])
        nul.append(cross(XA, MA, XB[:, p], MB[:, p]))
obs, nul = np.concatenate(obs), np.concatenate(nul)
n_ind = data[0][1].shape[1]
print(f"{os.path.basename(samples)} n={n_ind} blocks={len(data)} pairs={obs.size}")
print(f"  observed unlinked mean r2 {obs.mean():.4f}  median {np.median(obs):.4f}  share>0.1 {np.mean(obs > 0.1):.4f}")
print(f"  permuted unlinked mean r2 {nul.mean():.4f}  median {np.median(nul):.4f}  share>0.1 {np.mean(nul > 0.1):.4f}")
