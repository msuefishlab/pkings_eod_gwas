#!/usr/bin/env python3
"""ext_panel_chrom.py -- observed panel, one chromosome: panel_lib.segment_stats per mask block, and the
interchromosomal SNP draw (panel_lib.draw_interchrom) per sample.

Usage: ext_panel_chrom.py <blocks_npz_dir> <mask_blocks.tsv> <samples_def.tsv> <pairs_def.tsv> <chrom> <out_dir>
  samples_def.tsv: name <tab> comma-separated sample IDs;  pairs_def.tsv: bp_name <tab> tp_name
Writes <out_dir>/seg.<chrom>.tsv (chrom, block, kind, name, bin, value), <out_dir>/inter.<chrom>.npz,
<out_dir>/check.<chrom>.tsv (genotype vs site-table AN/AC agreement for the A1w samples).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_lib as P

bdir, mask_fn, sdef, pdef, chrom, out = sys.argv[1:7]
sets = {}
for line in open(sdef):
    if line.strip():
        n, ids = line.rstrip("\n").split("\t")
        sets[n] = ids.split(",")
pairs = [tuple(l.rstrip("\n").split("\t")) for l in open(pdef) if l.strip()]
blocks = []
for line in open(mask_fn):
    f = line.rstrip("\n").split("\t")
    if f[0] == chrom:
        blocks.append(int(f[1]))
rows = {n: ([], []) for n in sets}
os.makedirs(out, exist_ok=True)
tmp = os.path.join(out, f"seg.{chrom}.tsv.tmp")
chk = {"in_table": 0, "an_b_match": 0, "ac_b_match": 0, "an_t_match": 0, "ac_t_match": 0}
with open(tmp, "w") as fo:
    for b in blocks:
        z = np.load(os.path.join(bdir, f"{chrom}_{b}.npz"))
        names = list(z["samples"])
        col = {s: i for i, s in enumerate(names)}
        idx = {n: np.array([col[s] for s in ids]) for n, ids in sets.items()}
        G = z["G"]
        if G.size == 0:
            G = np.zeros((0, len(names)), np.int8)
        st = P.segment_stats(G, z["pos"], float(z["f"]), int(z["L"]), z["callable"], z["pol"], idx, pairs)
        for (kind, name), v in st.items():
            if np.ndim(v):
                for i, x in enumerate(v):
                    fo.write(f"{chrom}\t{b}\t{kind}\t{name}\t{i}\t{x:.10g}\n")
            else:
                fo.write(f"{chrom}\t{b}\t{kind}\t{name}\t\t{v:.10g}\n")
        for n, ix in idx.items():
            keep, X, M = P.ld_keep(G[:, ix])
            rows[n][0].append(X[keep]); rows[n][1].append(M[keep])
        if "BAM" in idx and "TP1" in idx and G.shape[0]:
            h = z["in_table"]
            for tag, ix, an, ac in (("b", idx["BAM"], z["an_b"], z["ac_b"]), ("t", idx["TP1"], z["an_t"], z["ac_t"])):
                Gs = G[h][:, ix]
                n_ = 2 * (Gs >= 0).sum(1); a_ = np.where(Gs >= 0, Gs, 0).sum(1)
                chk[f"an_{tag}_match"] += int(np.sum(n_ == an[h])); chk[f"ac_{tag}_match"] += int(np.sum(a_ == ac[h]))
            chk["in_table"] += int(h.sum())
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
with open(os.path.join(out, f"check.{chrom}.tsv"), "w") as fo:
    fo.write("\t".join(chk) + "\n" + "\t".join(str(v) for v in chk.values()) + "\n")
print(chrom, len(blocks), "blocks", chk, flush=True)
