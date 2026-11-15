#!/usr/bin/env python3
"""Per-file scan: count complete-data SNPs with BP-biased delta >= the Peak 4 / Peak 5 observed values
in non-overlapping windows of fixed size, for the A1w (Bambomo vs TP1) contrast. Peak regions kept."""
import gzip, sys, csv, collections
FN = sys.argv[1]; OUT = sys.argv[2]
THR = {"peak4": 124/134 - 11/94, "peak5": 105/134 - 0.0}
WINS = {"w354k": 353902, "w287k": 287163}
EPS = 1e-9
cnt = collections.Counter(); den = collections.Counter()
hdr = None
with gzip.open(FN, "rt") as fh:
    for r in csv.reader(fh, delimiter="\t"):
        if r[0].startswith("#"):
            hdr = {k.lstrip("#"): i for i, k in enumerate(r)}; continue
        an_b = int(r[hdr["AN_wBAM"]]); an_t = int(r[hdr["AN_wTP1"]])
        if an_b != 134 or an_t != 94:
            continue
        aa = r[hdr["aa"]]
        if aa not in ("R", "A"):
            continue
        ac_b = int(r[hdr["AC_wBAM"]]); ac_t = int(r[hdr["AC_wTP1"]])
        d_b = ac_b if aa == "R" else an_b - ac_b
        d_t = ac_t if aa == "R" else an_t - ac_t
        delta = d_b / an_b - d_t / an_t
        ch = r[0]; pos = int(r[1])
        for wn, w in WINS.items():
            key = (wn, ch, pos // w)
            den[key] += 1
            for t, thr in THR.items():
                if delta >= thr - EPS:
                    cnt[(wn, ch, pos // w, t)] += 1
with open(OUT, "w") as fh:
    fh.write("window_set\tchrom\twindow_index\tn_sites\tn_peak4\tn_peak5\n")
    for (wn, ch, wi), nd in sorted(den.items()):
        fh.write(f"{wn}\t{ch}\t{wi}\t{nd}\t{cnt.get((wn,ch,wi,'peak4'),0)}\t{cnt.get((wn,ch,wi,'peak5'),0)}\n")
