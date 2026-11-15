#!/usr/bin/env python3
"""08_peaks.py -- is Bambomo ancestry enriched at the six BP GWAS peaks?

Posteriors: hap_primary/post_all.npz (all markers, incl. peaks and the chr4 inversion; each fish at its own
genome-wide ML T and m). For a group of fish, the statistic in an interval is the mean over the AIMs inside it of
the group-mean posterior P(BAM) over the group's haplotypes. Peak interval = extended peak (peaks.bed cols 2-3).
Null: non-overlapping windows of the peak's length tiled from each chromosome start, wholly outside every peak
+/- 250 kb and the inversion, with >= 1 AIM. p_enrich = (1 + null >= peak)/(1 + n), one-sided; p_deplete the
same downward. Enriched: p <= 0.0083 (0.05/6); nominal: <= 0.05. Tested group: the ten 2019 Apassa wobbles.
Reported, not tested: the same for the control groups, and the number of haplotypes with mean posterior >= 0.5
in the peak.

Writes <out_dir>/results/peaks.tsv.

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/08_peaks.py <out_dir>
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import la_lib as L

GROUPS = {
    "APA2019_Wobble": lambda g, r: g == "APA2019_Wobble",
    "APA2009_TP_holdout": lambda g, r: g == "APA2009_TP" and "holdout_APA" in r,
    "APA2009_TP_q>=0.05": lambda g, r: g == "APA2009_TP" and "holdout_APA" not in r,
    "APA_6681": None, "APA_6713": None,
    "BEN2023_TP_holdout": lambda g, r: g == "BEN2023_TP",
    "BAM2009_holdout": lambda g, r: g.startswith("BAM2009"),
    "BAM2019_Wobble": lambda g, r: g == "BAM2019_Wobble",
}


def main():
    out = sys.argv[1]
    info = {}
    for line in open(f"{out}/samples.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "sample":
            info[f[0]] = (f[2], set(f[5].split(",")))
    z = np.load(f"{out}/hap_primary/post_all.npz")
    post, chrom, pos, hs = z["post"].astype(float), z["chrom"], z["pos"], z["hap_sample"]
    excl = L.excluded_regions()
    peaks = []
    for line in open("input_data/10_Introgression/peaks.bed"):
        f = line.split()
        peaks.append((int(f[3]), int(f[0]), int(f[1]), int(f[2])))
    peaks.sort()
    fai = {l.split()[0]: int(l.split()[1]) for l in open(
        "input_data/00_Reference_Genome/jordan-mic4273-mb-hirise-20xx3__06-13-2024__final_assembly.fasta.renamed.fa.fai")}
    rows = ["group\tn_haps\tpeak\tchrom\tstart\tend\tn_aims\tpeak_stat\tn_haps_post_ge_0.5\tnull_n\tnull_median\t"
            "null_p95\tp_enrich\tp_deplete\tcall"]
    for gname, sel in GROUPS.items():
        if sel is None:
            cols = np.flatnonzero(hs == gname)
        else:
            cols = np.flatnonzero([sel(*info[s]) for s in hs])
        gm = post[:, cols].mean(1)
        for pid, c, a, b in peaks:
            inside = (chrom == c) & (pos >= a) & (pos <= b)
            Lk = b - a
            if inside.sum() == 0:
                rows.append(f"{gname}\t{cols.size}\t{pid}\tchr{c}\t{a}\t{b}\t0\tNA\tNA\t0\tNA\tNA\tNA\tNA\tnot_testable")
                continue
            stat = gm[inside].mean()
            nh = int((post[inside][:, cols].mean(0) >= 0.5).sum())
            null = []
            for cc in range(1, 26):
                ch = f"chr{cc}"
                on = chrom == cc
                pc, gc = pos[on], gm[on]
                for w0 in range(0, fai[ch] - Lk + 1, Lk):
                    w1 = w0 + Lk
                    if any(x == ch and w0 < hi and w1 > lo for x, lo, hi in excl):
                        continue
                    i0, i1 = np.searchsorted(pc, w0), np.searchsorted(pc, w1)
                    if i1 > i0:
                        null.append(gc[i0:i1].mean())
            null = np.array(null)
            pe = (1 + (null >= stat).sum()) / (1 + null.size)
            pd = (1 + (null <= stat).sum()) / (1 + null.size)
            call = "enriched" if pe <= 0.05 / 6 else ("nominal" if pe <= 0.05 else "not_enriched")
            rows.append(f"{gname}\t{cols.size}\t{pid}\tchr{c}\t{a}\t{b}\t{inside.sum()}\t{stat:.4f}\t{nh}\t{null.size}\t"
                        f"{np.median(null):.4f}\t{np.percentile(null, 95):.4f}\t{pe:.4g}\t{pd:.4g}\t"
                        f"{call if gname == 'APA2019_Wobble' else 'descriptive'}")
    os.makedirs(f"{out}/results", exist_ok=True)
    open(f"{out}/results/peaks.tsv", "w").write("\n".join(rows) + "\n")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
