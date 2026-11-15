#!/usr/bin/env python3
"""11_apassa_pca.py -- within-Apassa PCA: where do the 2019-018/021 wobbles sit? (post hoc)

Genotypes: the unphased genotype tables written by simulation/26a_extract_temporal.sh (26_temporal/geno; PASS biallelic SNPs, FMT/DP >= 8 and
GQ >= 20, 2-kb regions every 20 kb, peaks +/- 250 kb and the chr4 inversion excluded). The 30 Apassa fish only
(metadata POP APA). SNPs with Apassa MAF >= 0.05 and no missing genotype among the 30; dosages centred and
scaled per SNP; PCA by SVD. Per fish: PC1-PC4 and the per-fish call rate over all extracted Apassa records (a
depth proxy, since DP < 8 calls are missing). Per PC: Pearson r with the 26c Bambomo score and with call rate,
and an exact two-sided Mann-Whitney test of 2019 wobbles against 2009 TP. For PCs with p <= 0.01, whether each
2019-018/021 fish and each 2009 wobble is nearer the 2019-017 centroid or the 2009 TP centroid.

Writes <out_dir>/results/apassa_pca.tsv and results/apassa_pca_summary.txt.

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/11_apassa_pca.py <out_dir>
"""
import gzip, sys
import numpy as np
from scipy.stats import mannwhitneyu, pearsonr

D = "output_data/11_Drift_vs_Selection/simulation/26_temporal"
CODE = {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2, "0|0": 0, "0|1": 1, "1|0": 1, "1|1": 2}


def main():
    O = sys.argv[1]
    names = open(f"{D}/samples.txt").read().split()
    info = {}
    for line in open(f"{O}/samples.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "sample":
            info[f[0]] = (f[2], f[4], float(f[3]))
    cols = [i for i, s in enumerate(names) if s in info and info[s][0].startswith("APA")]
    apa = [names[i] for i in cols]
    assert len(apa) == 30, len(apa)
    rows, called, total = [], np.zeros(30), 0
    for c in range(1, 26):
        with gzip.open(f"{D}/geno/chr{c}.tsv.gz", "rt") as fh:
            for line in fh:
                g = line.rstrip("\n").split("\t")[2:]
                x = np.array([CODE.get(g[i], -1) for i in cols])
                ok = x >= 0
                called += ok; total += 1
                if not ok.all():
                    continue
                p = x.mean() / 2
                if 0.05 <= p <= 0.95 and x.std() > 0:          # drops all-heterozygote sites
                    rows.append(x)
    X = np.array(rows, float)
    X = (X - X.mean(1, keepdims=True)) / X.std(1, keepdims=True)
    U, S, Vt = np.linalg.svd(X.T, full_matrices=False)
    pcs = U[:, :4] * S[:4]
    ve = S ** 2 / (S ** 2).sum()
    cr = called / total
    grp = [info[s][0] for s in apa]; site = [info[s][1] for s in apa]; q = np.array([info[s][2] for s in apa])
    out = ["sample\tgroup\tsite\tq26c\tcall_rate\tPC1\tPC2\tPC3\tPC4"]
    for i, s in enumerate(apa):
        out.append(f"{s}\t{grp[i]}\t{site[i]}\t{q[i]:.3f}\t{cr[i]:.4f}\t" + "\t".join(f"{v:.3f}" for v in pcs[i]))
    open(f"{O}/results/apassa_pca.tsv", "w").write("\n".join(out) + "\n")
    w19 = np.array([g == "APA2019_Wobble" for g in grp]); tp = np.array([g == "APA2009_TP" for g in grp])
    s17 = np.array([s == "2019-017" for s in site])
    lines = [f"SNPs {X.shape[0]} (complete in 30 Apassa fish, MAF >= 0.05) of {total} extracted records"]
    for k in range(4):
        v = pcs[:, k]
        p = mannwhitneyu(v[w19], v[tp], alternative="two-sided", method="exact").pvalue
        lines.append(f"PC{k + 1}: {100 * ve[k]:.2f}% var; r with q26c {pearsonr(v, q)[0]:+.3f}; r with call rate "
                     f"{pearsonr(v, cr)[0]:+.3f}; 2019 wobbles vs 2009 TP Mann-Whitney p {p:.2g}")
        if p <= 0.01:
            c17, ctp = v[s17].mean(), v[tp].mean()
            for i, s in enumerate(apa):
                if (grp[i] == "APA2019_Wobble" and not s17[i]) or grp[i] == "APA2009_Wobble":
                    near = "2019-017" if abs(v[i] - c17) < abs(v[i] - ctp) else "2009_TP"
                    lines.append(f"    {s} ({grp[i]}, {site[i]}): PC{k + 1} {v[i]:+.2f}; 2019-017 centroid {c17:+.2f}, "
                                 f"2009 TP centroid {ctp:+.2f} -> nearer {near}")
    open(f"{O}/results/apassa_pca_summary.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(out)); print("\n".join(lines))


if __name__ == "__main__":
    main()
