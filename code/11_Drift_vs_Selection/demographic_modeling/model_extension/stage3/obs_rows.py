#!/usr/bin/env python3
"""obs_rows.py -- observed panel rows for the Stage 3 sensitivity treatments, one chromosome (format of ld_panel/
ext_panel_chrom.py, so ext_aggregate.py runs on the output unchanged).

Rows:
  primary   Stage 1 block genotypes unchanged (re-check: must reproduce the Stage 1 panel seg files exactly)
  callable  only SNVs at A1w-callable positions (Stage 1 site-table flag), same blocks, filters, bootstrap
  d019      genotypes set missing where DP < 8 or GQ < 20 (records_full.py bit 4), all SNVs
Statistics: panel_lib.segment_stats and panel_lib.draw_interchrom, unchanged, for BAM and TP1 (pair BAM|TP1).
Usage: obs_rows.py <row> <blocks_npz_dir> <rec_dir> <mask_blocks.tsv> <samples_def.tsv> <pairs_def.tsv> <chrom> <out_dir>
"""
import os, sys
import numpy as np
sys.path.insert(0, os.environ["PANEL_DIR"])
import panel_lib as P  # noqa: E402

row, bdir, rdir, mask_fn, sdef, pdef, chrom, out = sys.argv[1:9]
sets = {}
for line in open(sdef):
    if line.strip():
        n, ids = line.rstrip("\n").split("\t")
        sets[n] = ids.split(",")
pairs = [tuple(l.rstrip("\n").split("\t")) for l in open(pdef) if l.strip()]
blocks = [int(l.split("\t")[1]) for l in open(mask_fn) if l.split("\t")[0] == chrom]
rows = {n: ([], []) for n in sets}
os.makedirs(out, exist_ok=True)
tmp = os.path.join(out, f"seg.{chrom}.tsv.tmp")
with open(tmp, "w") as fo:
    for b in blocks:
        z = np.load(os.path.join(bdir, f"{chrom}_{b}.npz"))
        names = list(z["samples"])
        col = {s: i for i, s in enumerate(names)}
        idx = {n: np.array([col[s] for s in ids]) for n, ids in sets.items()}
        G = z["G"]
        if G.size == 0:
            G = np.zeros((0, len(names)), np.int8)
        pos, cal, pol = z["pos"], z["callable"], z["pol"]
        if row == "callable":
            G, pos, cal, pol = G[cal], pos[cal], cal[cal], pol[cal]
        elif row == "d019" and len(pos):
            r = np.load(os.path.join(rdir, f"{chrom}_{b}.rec.npz"))
            ru = list(r["samples"]); rc = np.array([ru.index(s) for s in names])     # records cols -> npz order
            j = np.searchsorted(r["rpos"], pos)
            assert np.all(r["rpos"][np.clip(j, 0, len(r["rpos"]) - 1)] == pos)
            d019 = ((r["code"][j][:, rc] >> 4) & 1).astype(bool)
            G = np.where(d019, -1, G).astype(np.int8)
        elif row != "primary" and row != "d019":
            sys.exit(f"unknown row {row}")
        st = P.segment_stats(G, pos, float(z["f"]), int(z["L"]), cal, pol, idx, pairs)
        for (kind, name), v in st.items():
            if np.ndim(v):
                for i, x in enumerate(v):
                    fo.write(f"{chrom}\t{b}\t{kind}\t{name}\t{i}\t{x:.10g}\n")
            else:
                fo.write(f"{chrom}\t{b}\t{kind}\t{name}\t\t{v:.10g}\n")
        for n, ix in idx.items():
            keep, X, M = P.ld_keep(G[:, ix])
            rows[n][0].append(X[keep]); rows[n][1].append(M[keep])
os.replace(tmp, os.path.join(out, f"seg.{chrom}.tsv"))
c = int(chrom[3:])
arr = {}
for n, (xs, ms) in rows.items():
    X = np.concatenate(xs) if xs else np.zeros((0, len(sets[n])), np.int8)
    M = np.concatenate(ms) if ms else np.zeros((0, len(sets[n])), np.int8)
    d = P.draw_interchrom({c: (X, M)})[c]
    arr[f"{n}__X"], arr[f"{n}__M"], arr[f"{n}__Xp"], arr[f"{n}__Mp"] = (a.astype(np.int8) for a in d)
    arr[f"{n}__m"] = np.array(X.shape[0])
np.savez_compressed(os.path.join(out, f"inter.{chrom}.tmp.npz"), **arr)
os.replace(os.path.join(out, f"inter.{chrom}.tmp.npz"), os.path.join(out, f"inter.{chrom}.npz"))
print(row, chrom, len(blocks), "blocks", flush=True)
