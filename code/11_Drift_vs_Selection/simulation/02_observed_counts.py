#!/usr/bin/env python3
"""Per-file scan of the drift-test site table (demographic_modeling/fixed_differences/): count SNPs with BP-biased delta >= the Peak 4 / Peak 5
observed values, genome-wide and inside the peak regions (whole-population Bambomo vs TP1)."""
import gzip, sys, csv, os
FN = sys.argv[1]; OUT = sys.argv[2]
REG = {"peak1": ("chr6", 29448847, 29509385), "peak4": ("chr16", 10728257, 11082159),
       "peak5": ("chr17", 11494062, 11781225)}
THR = {"peak4": 124/134 - 11/94, "peak5": 105/134 - 0.0}
EPS = 1e-9
res = {(pk, t): 0 for pk in REG for t in THR}
gw = {t: 0 for t in THR}; gw_out = {t: 0 for t in THR}
n_callable = 0; n_both_called = 0
hdr = None
with gzip.open(FN, "rt") as fh:
    for r in csv.reader(fh, delimiter="\t"):
        if r[0].startswith("#"):
            hdr = {k.lstrip("#"): i for i, k in enumerate(r)}; continue
        if r[hdr["excl"]] != "0":
            continue
        n_callable += 1
        an_b = int(r[hdr["AN_wBAM"]]); ac_b = int(r[hdr["AC_wBAM"]])
        an_t = int(r[hdr["AN_wTP1"]]); ac_t = int(r[hdr["AC_wTP1"]])
        if an_b == 0 or an_t == 0:
            continue
        n_both_called += 1
        aa = r[hdr["aa"]]
        d_b = ac_b if aa == "R" else an_b - ac_b
        d_t = ac_t if aa == "R" else an_t - ac_t
        delta = d_b / an_b - d_t / an_t
        ch = r[0]; pos = int(r[1]); inreg = None
        for pk, (c, s, e) in REG.items():
            if ch == c and s <= pos <= e:
                inreg = pk; break
        for t, thr in THR.items():
            if delta >= thr - EPS:
                gw[t] += 1
                if inreg: res[(inreg, t)] += 1
                else: gw_out[t] += 1
with open(OUT, "w") as fh:
    fh.write("scope\tthreshold\tn\n")
    fh.write(f"callable_snps\tNA\t{n_callable}\n")
    fh.write(f"both_called_snps\tNA\t{n_both_called}\n")
    for t in THR:
        fh.write(f"genome_all\t{t}\t{gw[t]}\n")
        fh.write(f"genome_outside_peaks\t{t}\t{gw_out[t]}\n")
        for pk in REG:
            fh.write(f"region_{pk}\t{t}\t{res[(pk,t)]}\n")
