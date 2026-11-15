#!/usr/bin/env python3
"""ext_sitetable.py -- per-SNP callable flag and polarization from one drift-test site-table chunk (A1w rule).

callable (A1w, exactly the 03_spectra.py --mode spectra rule behind L_scan of sfs/A1w.block_scan.tsv):
  ExcHet >= 0.001 in every whole population (paralog filter), AN_wBAM >= n_bp and AN_wTP1 >= n_tp
  (sfs/projections.tsv, A1w: 60 / 60), polarized (aa R or A), and the 100-kb bin not an A1w call desert
  (sfs/deserts.tsv).
pol: +1 if aa == R (derived = ALT), -1 if aa == A (derived = REF), 0 otherwise.
Also keeps AN/AC of wBAM and wTP1 for the genotype cross-check.

Usage: ext_sitetable.py <snps.N.tsv.gz> <projections.tsv> <deserts.tsv> <out.npz>
"""
import gzip, sys
import numpy as np

fn, proj_fn, des_fn, out = sys.argv[1:5]
n1 = n2 = None
for line in open(proj_fn):
    f = line.rstrip("\n").split("\t")
    if f[0] == "A1w":
        n1, n2 = int(f[1]), int(f[2])
des = set()
for line in open(des_fn):
    f = line.rstrip("\n").split("\t")
    if f[0] == "A1w":
        des.add((f[1], int(f[2])))
to_int = lambda x: 0 if x in (".", "") else int(x)
chrom, pos, cal, pol, anb, acb, ant, act = [], [], [], [], [], [], [], []
with gzip.open(fn, "rt") as fh:
    hdr = fh.readline().lstrip("#").rstrip("\n").split("\t")
    ix = {c: i for i, c in enumerate(hdr)}
    exc = [ix[c] for c in hdr if c.startswith("ExcHet_")]
    ib, jb, it, jt = ix["AN_wBAM"], ix["AC_wBAM"], ix["AN_wTP1"], ix["AC_wTP1"]
    for line in fh:
        f = line.rstrip("\n").split("\t")
        aa = f[4]
        p = 1 if aa == "R" else (-1 if aa == "A" else 0)
        ok_exc = not any(f[i] not in (".", "") and float(f[i]) < 0.001 for i in exc)
        Nb, Nt = to_int(f[ib]), to_int(f[it])
        x = int(f[1]) - 1
        c = ok_exc and p != 0 and Nb >= n1 and Nt >= n2 and (f[0], x // 100_000) not in des
        chrom.append(f[0]); pos.append(int(f[1])); cal.append(c); pol.append(p)
        anb.append(Nb); acb.append(to_int(f[jb])); ant.append(Nt); act.append(to_int(f[jt]))
tmp = out[:-4] + ".tmp.npz"
np.savez_compressed(tmp, chrom=np.array(chrom), pos=np.array(pos, np.int64), callable=np.array(cal, bool),
                    pol=np.array(pol, np.int8), an_b=np.array(anb, np.int16), ac_b=np.array(acb, np.int16),
                    an_t=np.array(ant, np.int16), ac_t=np.array(act, np.int16))
import os
os.replace(tmp, out)
print(fn, len(pos), "rows", int(np.sum(cal)), "callable")
