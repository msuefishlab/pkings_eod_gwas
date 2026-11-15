#!/usr/bin/env python3
"""24a_site_scan.py -- per-chunk extract for the empirical clustering test (24b).

Same site filter as 03_window_scan.py: complete data in the A1w samples (AN_wBAM = 134, AN_wTP1 = 94) and
polarized (aa R or A). Delta = BP-derived allele frequency in Bambomo minus that in TP1. Writes
  <out>.qual.tsv   chrom, pos, delta for every site with delta >= 0.78 (covers the Peak 4 and Peak 5 thresholds)
  <out>.bins.tsv   chrom, 1-kb bin, number of sites in the filtered set (SNP density)

Usage: 24a_site_scan.py <snps.N.tsv.gz> <out_prefix>
"""
import collections, csv, gzip, sys

FN, OUT = sys.argv[1], sys.argv[2]
DMIN = 0.78 - 1e-9
bins = collections.Counter()
with gzip.open(FN, "rt") as fh, open(OUT + ".qual.tsv", "w") as q:
    hdr = None
    for r in csv.reader(fh, delimiter="\t"):
        if r[0].startswith("#"):
            hdr = {k.lstrip("#"): i for i, k in enumerate(r)}
            continue
        an_b, an_t = int(r[hdr["AN_wBAM"]]), int(r[hdr["AN_wTP1"]])
        if an_b != 134 or an_t != 94:
            continue
        aa = r[hdr["aa"]]
        if aa not in ("R", "A"):
            continue
        ac_b, ac_t = int(r[hdr["AC_wBAM"]]), int(r[hdr["AC_wTP1"]])
        d_b = ac_b if aa == "R" else an_b - ac_b
        d_t = ac_t if aa == "R" else an_t - ac_t
        delta = d_b / an_b - d_t / an_t
        pos = int(r[1])
        bins[(r[0], pos // 1000)] += 1
        if delta >= DMIN:
            q.write(f"{r[0]}\t{pos}\t{delta:.6f}\n")
with open(OUT + ".bins.tsv", "w") as fh:
    for (ch, b), n in sorted(bins.items()):
        fh.write(f"{ch}\t{b}\t{n}\n")
