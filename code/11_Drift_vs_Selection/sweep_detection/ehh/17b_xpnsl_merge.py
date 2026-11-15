#!/usr/bin/env python3
"""
17b_xpnsl_merge.py
Merge the per-chromosome raw XP-nSL scores of 17_xpnsl.py and standardize them genome-wide per
contrast (mean 0, SD 1 over all 25 chromosomes, as selscan's norm).

Writes ehh/xpnsl/<BP>.genome_wide.txt.gz (chrom, pos, value), the layout load_derived() in
evidence/02_region_stats.R reads.

Usage: python3 17b_xpnsl_merge.py <raw_dir> <out_dir>
"""

import glob
import gzip
import os
import sys

import numpy as np

raw_dir, out_dir = sys.argv[1], sys.argv[2]
os.makedirs(out_dir, exist_ok=True)
for pair in ("BP1_TP1", "BP2_TP2", "BP3_TP2"):
    files = sorted(glob.glob(os.path.join(raw_dir, f"{pair}.chr*.xpnsl_raw.tsv.gz")))
    if len(files) != 25:
        sys.exit(f"{pair}: expected 25 chromosome files, found {len(files)}")
    chrom, pos, raw = [], [], []
    for f in files:
        with gzip.open(f, "rt") as fh:
            next(fh)
            for line in fh:
                c, p, v = line.split()
                chrom.append(c)
                pos.append(int(p))
                raw.append(float(v))
    raw = np.array(raw)
    mu, sd = raw.mean(), raw.std()
    z = (raw - mu) / sd
    bp = pair.split("_")[0]
    out = os.path.join(out_dir, f"{bp}.genome_wide.txt.gz")
    order = sorted(range(len(pos)), key=lambda i: (int(chrom[i][3:]), pos[i]))
    with gzip.open(out, "wt") as fh:
        fh.write("chrom\tpos\tvalue\n")
        for i in order:
            fh.write(f"{chrom[i]}\t{pos[i]}\t{z[i]:.6f}\n")
    print(f"{pair}: {len(raw)} sites, raw mean {mu:.4f}, sd {sd:.4f} -> {out}")
