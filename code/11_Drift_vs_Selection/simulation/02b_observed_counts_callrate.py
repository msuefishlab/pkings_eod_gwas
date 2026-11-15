#!/usr/bin/env python3
"""Per-file scan of the drift-test site table (demographic_modeling/fixed_differences/) for the A1w (whole-population Bambomo vs TP1) contrast.
Counts SNPs with BP-biased delta >= the observed Peak 4 / Peak 5 values, under two call-rate filters,
split by whether the site lies inside a GWAS peak region (excl flag / dbin) or not."""
import gzip, sys, csv
FN = sys.argv[1]; OUT = sys.argv[2]
REG = {"peak1": ("chr6", 29448847, 29509385), "peak4": ("chr16", 10728257, 11082159),
       "peak5": ("chr17", 11494062, 11781225)}
THR = {"peak4": 124/134 - 11/94, "peak5": 105/134 - 0.0}
FILTERS = {"complete": (134, 94), "call90": (120, 84)}
EPS = 1e-9
cnt = {}
def add(k):
    cnt[k] = cnt.get(k, 0) + 1
hdr = None
with gzip.open(FN, "rt") as fh:
    for r in csv.reader(fh, delimiter="\t"):
        if r[0].startswith("#"):
            hdr = {k.lstrip("#"): i for i, k in enumerate(r)}; continue
        an_b = int(r[hdr["AN_wBAM"]]); ac_b = int(r[hdr["AC_wBAM"]])
        an_t = int(r[hdr["AN_wTP1"]]); ac_t = int(r[hdr["AC_wTP1"]])
        if an_b == 0 or an_t == 0:
            continue
        aa = r[hdr["aa"]]
        if aa not in ("R", "A"):
            continue
        d_b = ac_b if aa == "R" else an_b - ac_b
        d_t = ac_t if aa == "R" else an_t - ac_t
        delta = d_b / an_b - d_t / an_t
        ch = r[0]; pos = int(r[1])
        inreg = ""
        for pk, (c, s, e) in REG.items():
            if ch == c and s <= pos <= e:
                inreg = pk; break
        scope = f"region_{inreg}" if inreg else "outside_peak_regions"
        for fname, (mb, mt) in FILTERS.items():
            if an_b >= mb and an_t >= mt:
                add((scope, fname, "denominator"))
                for t, thr in THR.items():
                    if delta >= thr - EPS:
                        add((scope, fname, t))
with open(OUT, "w") as fh:
    fh.write("scope\tfilter\tstat\tn\n")
    for (a, b, cc), n in sorted(cnt.items()):
        fh.write(f"{a}\t{b}\t{cc}\t{n}\n")
