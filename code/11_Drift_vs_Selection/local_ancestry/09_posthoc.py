#!/usr/bin/env python3
"""09_posthoc.py -- post hoc descriptive comparisons following 07_match.py and 08_peaks.py.

1. Long-tract share: F_long / q_m (the share of a fish's Bambomo-like ancestry that sits in called tracts
   >= 2 Mb) and the longest called tract, for every simulated cell (m = 0.10-0.20) and every real Apassa fish
   outside the panels. The excess share subtracts the 2009 Apassa TP median (F_long and q_m) from the fish.
2. Per-haplotype mean posterior P(BAM) inside each extended peak for the Apassa fish outside the panels
   (frequency-based; not a haplotype identity test).

Writes <out_dir>/results/posthoc_long_share.tsv and results/peak_copies.tsv.

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/09_posthoc.py <out_dir>
"""
import csv, sys
import numpy as np


def main():
    O = sys.argv[1]
    S = list(csv.DictReader(open(f"{O}/sim/fish.tsv"), delimiter="\t"))
    R = list(csv.DictReader(open(f"{O}/results/per_fish.tsv"), delimiter="\t"))
    rows = ["source\tid\tT\tm\tsite\tq_m\tF_long\tshare\tshare_p10\tshare_p90\texcess_share\tlongest_Mb\t"
            "longest_p10\tlongest_p90\tfrac_longest_ge_20Mb"]
    for m in ("0.1", "0.15", "0.2"):
        for T in ("1", "2", "3", "5", "10", "20", "30", "50", "100", "200", "300", "500", "1000"):
            rr = [r for r in S if r["T"] == T and r["m"] == m]
            sh = np.array([float(r["F_long"]) / float(r["q_m"]) for r in rr])
            lg = np.array([float(r["longest"]) for r in rr]) / 1e6
            rows.append(f"sim\tcell\t{T}\t{m}\t-\t{np.median([float(r['q_m']) for r in rr]):.4f}\t"
                        f"{np.median([float(r['F_long']) for r in rr]):.4f}\t{np.median(sh):.3f}\t"
                        f"{np.percentile(sh, 10):.3f}\t{np.percentile(sh, 90):.3f}\tNA\t{np.median(lg):.2f}\t"
                        f"{np.percentile(lg, 10):.2f}\t{np.percentile(lg, 90):.2f}\t{np.mean(lg >= 20):.2f}")
    base = [r for r in R if r["group"] == "APA2009_TP"]
    bF = np.median([float(r["F_long"]) for r in base]); bq = np.median([float(r["q_m"]) for r in base])
    rows.append(f"baseline\tAPA2009_TP_median_n{len(base)}\tNA\tNA\t2009-018\t{bq:.4f}\t{bF:.4f}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")
    for r in R:
        if r["group"].startswith("APA"):
            q, F = float(r["q_m"]), float(r["F_long"])
            ex = (F - bF) / (q - bq) if q - bq > 0.02 else float("nan")
            rows.append(f"real\t{r['sample']}\tNA\tNA\t{r['site']}\t{q:.4f}\t{F:.4f}\t{F / q if q > 0 else float('nan'):.3f}"
                        f"\tNA\tNA\t{ex:.3f}\t{float(r['longest']) / 1e6:.2f}\tNA\tNA\tNA")
    open(f"{O}/results/posthoc_long_share.tsv", "w").write("\n".join(rows) + "\n")

    z = np.load(f"{O}/hap_primary/post_all.npz")
    post, chrom, pos, hs = z["post"], z["chrom"], z["pos"], z["hap_sample"]
    pk = {}
    for l in open("input_data/10_Introgression/peaks.bed"):
        f = l.split(); pk[int(f[3])] = (int(f[0]), int(f[1]), int(f[2]))
    out = ["sample\tgroup\tsite\t" + "\t".join(f"P{p}_h0\tP{p}_h1" for p in range(1, 7))]
    for r in R:
        if not r["group"].startswith("APA"):
            continue
        c = np.flatnonzero(hs == r["sample"])
        vals = []
        for p in range(1, 7):
            ch, a, b = pk[p]
            ins = (chrom == ch) & (pos >= a) & (pos <= b)
            v = post[ins][:, c].mean(0); vals += [f"{v[0]:.3f}", f"{v[1]:.3f}"]
        out.append(f"{r['sample']}\t{r['group']}\t{r['site']}\t" + "\t".join(vals))
    open(f"{O}/results/peak_copies.tsv", "w").write("\n".join(out) + "\n")
    print("\n".join(rows[-13:]))


if __name__ == "__main__":
    main()
