#!/usr/bin/env python3
"""12_peak_runs.py -- Bambomo-like runs through the six BP peaks, per Apassa haplotype (post hoc).

Descriptive only, no test.

Input: hap_primary/post_all.npz (posterior P(BAM) at every marker, peaks included; the called tracts in
tracts.tsv exclude peaks +/- 250 kb, so they cannot be used here). For every Apassa fish (metadata group APA*)
and haplotype whose mean posterior over the markers inside an extended peak (input_data/10_Introgression/peaks.bed)
is >= 0.5: the run of consecutive markers with P(BAM) >= 0.5 that contains the peak marker with the highest
posterior. Runs are frequency-based ancestry calls, not a haplotype identity test (similar is not the same).

Writes <out_dir>/results/peak_runs.tsv (one row per carrier haplotype) and results/peak_runs_summary.tsv (per
peak and group: fish with 0, 1 or 2 Bambomo-like copies).

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/12_peak_runs.py <out_dir>
"""
import sys
import numpy as np


def main():
    O = sys.argv[1]
    info = {}
    for line in open(f"{O}/samples.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "sample":
            info[f[0]] = (f[2], f[4], f[5])
    z = np.load(f"{O}/hap_primary/post_all.npz")
    post, chrom, pos, hs = z["post"], z["chrom"], z["pos"], z["hap_sample"]
    pk = {}
    for line in open("input_data/10_Introgression/peaks.bed"):
        f = line.split()
        pk[int(f[3])] = (int(f[0]), int(f[1]), int(f[2]))
    fish = [s for s in dict.fromkeys(hs) if s in info and info[s][0].startswith("APA") and info[s][2] == "-"
            or s in info and info[s][0].startswith("APA") and "holdout" in info[s][2]]
    rows = ["peak\tsample\tgroup\tsite\thap\tpeak_mean_post\trun_start\trun_end\trun_Mb\trun_markers"]
    copies = {}
    for p in sorted(pk):
        c, a, b = pk[p]
        on = np.flatnonzero(chrom == c)
        P = pos[on]
        ins = (P >= a) & (P <= b)
        for s in fish:
            n = 0
            for h, j in enumerate(np.flatnonzero(hs == s)):
                v = post[on, j]
                m = v[ins].mean()
                if m < 0.5:
                    continue
                n += 1
                i0 = np.flatnonzero(ins)[np.argmax(v[ins])]
                L = i0
                while L > 0 and v[L - 1] >= 0.5:
                    L -= 1
                R = i0
                while R < len(v) - 1 and v[R + 1] >= 0.5:
                    R += 1
                rows.append(f"{p}\t{s}\t{info[s][0]}\t{info[s][1]}\t{h}\t{m:.3f}\tchr{c}:{P[L]}\tchr{c}:{P[R]}\t"
                            f"{(P[R] - P[L]) / 1e6:.2f}\t{R - L + 1}")
            copies[(p, s)] = n
    open(f"{O}/results/peak_runs.tsv", "w").write("\n".join(rows) + "\n")
    groups = [("APA2019_Wobble 2019-017", lambda s: info[s][0] == "APA2019_Wobble" and info[s][1] == "2019-017"),
              ("APA2019_Wobble 2019-018/021", lambda s: info[s][0] == "APA2019_Wobble" and info[s][1] != "2019-017"),
              ("APA2009_Wobble", lambda s: info[s][0] == "APA2009_Wobble"),
              ("APA2009_TP outside panels", lambda s: info[s][0] == "APA2009_TP")]
    out = ["peak\tgroup\tn_fish\tfish_0_copies\tfish_1_copy\tfish_2_copies\tcarriers"]
    for p in sorted(pk):
        for g, sel in groups:
            fs = [s for s in fish if sel(s)]
            k = [copies[(p, s)] for s in fs]
            car = ",".join(f"{s}({copies[(p, s)]})" for s in fs if copies[(p, s)])
            out.append(f"{p}\t{g}\t{len(fs)}\t{k.count(0)}\t{k.count(1)}\t{k.count(2)}\t{car or '-'}")
    open(f"{O}/results/peak_runs_summary.tsv", "w").write("\n".join(out) + "\n")
    print("\n".join(rows)); print("\n".join(out))


if __name__ == "__main__":
    main()
