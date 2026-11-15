#!/usr/bin/env python3
"""genome_post.py -- one simulated (or observed-row) genome directory <d> with seg/ (25 seg.<chrom>.tsv + inter.<chrom>.npz):
runs ld_panel/ext_aggregate.py UNCHANGED (panel.tsv, ld_curves.tsv: genome / calibration / heldout point values and
chromosome-bootstrap bands), and writes chromld.npz = per-chromosome LD sums (n_pairs, sum_r2 per bin) per sample, so
r/mu can be refitted with chromosome-bootstrap intervals (23e) without re-reading the seg files.
Usage: genome_post.py <d> <chrom_split.tsv> <samples_def.tsv> <pairs_def.tsv>
"""
import os, subprocess, sys
import numpy as np

d, split, sdef, pdef = sys.argv[1:5]
seg = os.path.join(d, "seg")
if not (os.path.exists(os.path.join(d, "panel.tsv")) and os.path.exists(os.path.join(d, "ld_curves.tsv"))):
    subprocess.run([sys.executable, os.path.join(os.environ["PANEL_DIR"], "ext_aggregate.py"), seg, split, sdef, pdef,
                    os.path.join(d, "panel.tsv"), os.path.join(d, "ld_curves.tsv")], check=True,
                   stdout=subprocess.DEVNULL)
names = [l.split("\t")[0] for l in open(sdef) if l.strip()]
acc = {}
for fn in sorted(os.listdir(seg)):
    if not (fn.startswith("seg.") and fn.endswith(".tsv")):
        continue
    ch = fn[4:-4]
    A = {(n, k): np.zeros(50) for n in names for k in ("ld_cnt", "ld_tot")}
    for line in open(os.path.join(seg, fn)):
        c, b, kind, name, i, v = line.rstrip("\n").split("\t")
        if kind in ("ld_cnt", "ld_tot") and name in names:
            A[(name, kind)][int(i)] += float(v)
    for (n, k), v in A.items():
        acc[f"{ch}|{n}|{k}"] = v
tmp = os.path.join(d, "chromld.tmp.npz")
np.savez_compressed(tmp, **acc)
os.replace(tmp, os.path.join(d, "chromld.npz"))
print("post", d, flush=True)
