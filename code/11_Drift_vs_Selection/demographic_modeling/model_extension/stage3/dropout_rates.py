#!/usr/bin/env python3
"""dropout_rates.py -- forward heterozygote-dropout rate d(i, b) per individual i x observed-DP bin b.

  H_miss(i, b)   = sum over i's homozygous calls in bin b of 10^(-GQ/10)
  H_called(i, b) = number of i's heterozygous calls in bin b
  d(i, b)        = H_miss / (H_called + H_miss)            (0 when both are 0)
The sums run over the ld23-eligible SNPs of the 860 M-sum blocks: for a BAM fish the SNPs that pass ld23_lib.filter_snps
(MAF >= 0.10, call >= 0.90) in BAM, for a TP1 fish those passing in TP1 (panel_lib.ld_keep, unchanged), on the Stage 1
block genotypes (23a coding). GQ and the DP bin come from records_full.py at the same position. Homozygous calls with
GQ '.' contribute nothing and are counted (n_hom_gq_missing).
Usage: dropout_rates.py <blocks_npz_dir> <rec_dir> <mask_blocks.tsv> <samples_union.txt> <samples_def.tsv> <out.tsv> <nproc>
"""
import csv, multiprocessing as mp, os, sys
import numpy as np

bdir, rdir, mask_fn, union_fn, sdef, out, nproc = sys.argv[1:8]
sys.path.insert(0, os.environ["PANEL_DIR"])
import panel_lib as P  # noqa: E402

union = [l.strip() for l in open(union_fn) if l.strip()]
sets = {}
for line in open(sdef):
    if line.strip():
        nm, ids = line.rstrip("\n").split("\t")
        if nm in ("BAM", "TP1"):
            sets[nm] = ids.split(",")
SIDX = {nm: np.array([union.index(x) for x in ids]) for nm, ids in sets.items()}
blocks = [(r["chrom"], int(r["block"])) for r in csv.DictReader(open(mask_fn), delimiter="\t")]
NB = 5


def one(cb):
    ch, b = cb
    z = np.load(os.path.join(bdir, f"{ch}_{b}.npz"))
    zn = list(z["samples"]); zc = np.array([zn.index(x) for x in union])
    G = z["G"][:, zc] if z["G"].size else np.zeros((0, len(union)), np.int8)
    zpos = z["pos"]
    r = np.load(os.path.join(rdir, f"{ch}_{b}.rec.npz"))
    rpos, code, vrow, gq = r["rpos"], r["code"], r["vrow"], r["gq"]
    Hm = np.zeros((len(union), NB)); Hc = np.zeros((len(union), NB)); nh = np.zeros((len(union), NB))
    nmiss = np.zeros(len(union))
    if len(zpos) == 0:
        return Hm, Hc, nh, nmiss, 0
    j = np.searchsorted(rpos, zpos)
    assert np.all(rpos[np.clip(j, 0, len(rpos) - 1)] == zpos), f"{ch}:{b} npz SNV not in records"
    k = np.searchsorted(vrow, j)
    assert np.all(vrow[np.clip(k, 0, len(vrow) - 1)] == j), f"{ch}:{b} npz SNV without GQ row"
    nsnp = 0
    for nm, idx in SIDX.items():
        keep, _, _ = P.ld_keep(G[:, idx])
        if not keep.any():
            continue
        nsnp += int(keep.sum())
        Gs = G[keep][:, idx]
        bins = (code[j[keep]][:, idx] & 7).astype(int)
        q = gq[k[keep]][:, idx].astype(float)
        hom = (Gs == 0) | (Gs == 2)
        het = Gs == 1
        qok = q < 255
        w = np.where(hom & qok, 10.0 ** (-q / 10.0), 0.0)
        for bb in range(NB):
            m = bins == bb
            Hm[idx, bb] += (w * m).sum(0)
            Hc[idx, bb] += (het & m).sum(0)
            nh[idx, bb] += (hom & m).sum(0)
        nmiss[idx] += (hom & ~qok).sum(0)
    return Hm, Hc, nh, nmiss, nsnp


with mp.Pool(int(nproc)) as pool:
    res = pool.map(one, blocks, chunksize=4)
Hm = sum(r[0] for r in res); Hc = sum(r[1] for r in res); nh = sum(r[2] for r in res); nm = sum(r[3] for r in res)
nsnp = sum(r[4] for r in res)
LAB = ["0-3", "4-7", "8-11", "12-19", ">=20"]
memb = {i: nm_ for nm_, idx in SIDX.items() for i in idx}
with open(out + ".tmp", "w") as fo:
    fo.write("i\tsample_id\tsample\tbin\tdp_bin\tH_called\tH_miss\tn_hom\tn_hom_gq_missing\td\n")
    for i, sid in enumerate(union):
        for bb in range(NB):
            den = Hc[i, bb] + Hm[i, bb]
            d = Hm[i, bb] / den if den > 0 else 0.0
            fo.write(f"{i}\t{sid}\t{memb.get(i, '')}\t{bb}\t{LAB[bb]}\t{int(Hc[i, bb])}\t{Hm[i, bb]:.6g}\t{int(nh[i, bb])}\t"
                     f"{int(nm[i]) if bb == 0 else ''}\t{d:.8g}\n")
os.replace(out + ".tmp", out)
print("blocks", len(blocks), "ld23-eligible SNP x sample rows", nsnp, flush=True)
