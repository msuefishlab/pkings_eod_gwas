#!/usr/bin/env python3
"""prep_cache.py -- unpack records_full.py npz files into raw .npy (rpos, code, cal) in a node-local cache dir.
Usage: prep_cache.py <rec_dir> <mask_blocks.tsv> <cache_dir> <nproc>"""
import csv, multiprocessing as mp, os, sys
import numpy as np

rdir, mask_fn, cache, nproc = sys.argv[1:5]
os.makedirs(cache, exist_ok=True)
blocks = [(r["chrom"], int(r["block"])) for r in csv.DictReader(open(mask_fn), delimiter="\t")]


def one(cb):
    ch, b = cb
    p = os.path.join(cache, f"{ch}_{b}")
    if all(os.path.exists(p + x) for x in (".rpos.npy", ".code.npy", ".cal.npy")):
        return 0
    z = np.load(os.path.join(rdir, f"{ch}_{b}.rec.npz"))
    for k, x in (("rpos", ".rpos.npy"), ("code", ".code.npy"), ("callable", ".cal.npy")):
        np.save(p + ".tmp" + x, z[k]); os.replace(p + ".tmp" + x, p + x)
    return 1


with mp.Pool(int(nproc)) as pool:
    n = sum(pool.map(one, blocks, chunksize=2))
print("cache", cache, "blocks", len(blocks), "written", n, flush=True)
