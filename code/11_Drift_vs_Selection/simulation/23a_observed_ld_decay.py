#!/usr/bin/env python3
"""23a_observed_ld_decay.py -- observed LD decay for the r/mu fit (LD-decay test).

Segments are 200 of the drift test's 1-Mb scan blocks (sfs/A1w.block_scan.tsv) drawn at random (seed 20261030) from those with callable sequence
that lie wholly inside their chromosome and clear of the six extended peaks +/-250 kb and the chr4 inversion-like
polymorphism. Each block's callable fraction f = L_scan / 1 Mb converts physical distance to the simulation's
callable-layout distance. Genotypes come from the unphased VCF for one sample list (Bambomo = the 67 A1w fish; TP1
= the 47 APA + BEN fish), biallelic SNVs only; the statistic is ld23_lib.binned_r2.

Usage: 23a_observed_ld_decay.py <vcf> <samples.txt> <label> <block_scan.tsv> <peaks.bed> <out.tsv>
       [--blocks chrom:block,...]   score these blocks instead of the random draw (dense-window check)
Output: one row per segment x bin: label, chrom, block, f, n_snps, bin, n_pairs, sum_r2.
"""
import os, subprocess, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ld23_lib as L

N_SEG = 200
SEED = 20261030
FLANK = 250_000
CHR4_INVERSION = ("chr4", 40_880_000, 42_770_000)   # as fixed_differences/00_prepare_inputs.py
CHROMS = {f"chr{i}" for i in range(1, 26)}


def chrom_lengths(vcf):
    out = subprocess.run(["bcftools", "view", "-h", vcf], capture_output=True, text=True, check=True).stdout
    ln = {}
    for line in out.splitlines():
        if line.startswith("##contig=<ID="):
            body = line[len("##contig=<"):-1]
            kv = dict(x.split("=", 1) for x in body.split(","))
            ln[kv["ID"]] = int(kv["length"])
    return ln


def eligible_blocks(block_scan, peaks_bed, lengths):
    excl = [CHR4_INVERSION]
    for line in open(peaks_bed):
        p = line.split()
        excl.append((f"chr{p[0]}", int(p[1]) - FLANK, int(p[2]) + FLANK))
    out = []
    with open(block_scan) as fh:
        next(fh)
        for line in fh:
            ch, b, lscan = line.split("\t")[:3]
            b, lscan = int(b), int(lscan)
            s, e = b * 1_000_000, (b + 1) * 1_000_000
            if ch not in CHROMS or lscan <= 0 or e > lengths.get(ch, 0):
                continue
            if any(c == ch and s < xe and e > xs for c, xs, xe in excl):
                continue
            out.append((ch, b, lscan / 1_000_000))
    return out


def genotypes(vcf, samples, ch, s, e):
    cmd = ["bcftools", "view", "-r", f"{ch}:{s + 1}-{e}", "-S", samples, "-v", "snps", "-m2", "-M2", vcf]
    q = ["bcftools", "query", "-f", "%POS[\t%GT]\n"]
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    out = subprocess.run(q, stdin=p1.stdout, capture_output=True, text=True, check=True).stdout
    p1.wait()
    code = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/0": 1, "1/1": 2, "1|1": 2}
    pos, G = [], []
    for line in out.splitlines():
        f = line.split("\t")
        pos.append(int(f[0]))
        G.append([code.get(g, -1) for g in f[1:]])
    G = np.asarray(G, np.int8).reshape(len(pos), -1) if pos else np.empty((0, 0), np.int8)
    return np.asarray(pos, np.int64), G


def main():
    vcf, samples, label, block_scan, peaks_bed, out = sys.argv[1:7]
    chosen = None
    if "--blocks" in sys.argv:
        chosen = [(x.split(":")[0], int(x.split(":")[1])) for x in sys.argv[sys.argv.index("--blocks") + 1].split(",")]
    lengths = chrom_lengths(vcf)
    if chosen is None:
        el = eligible_blocks(block_scan, peaks_bed, lengths)
        rng = np.random.default_rng(SEED)
        segs = [el[i] for i in sorted(rng.choice(len(el), N_SEG, replace=False))]
    else:
        fr = {}
        with open(block_scan) as fh:
            next(fh)
            for line in fh:
                ch, b, lscan = line.split("\t")[:3]
                fr[(ch, int(b))] = int(lscan) / 1_000_000
        segs = [(ch, b, fr[(ch, b)]) for ch, b in chosen]
    with open(out, "w") as fo:
        fo.write("label\tchrom\tblock\tf\tn_snps\tbin\tn_pairs\tsum_r2\n")
        for k, (ch, b, f) in enumerate(segs, 1):
            pos, G = genotypes(vcf, samples, ch, b * 1_000_000, (b + 1) * 1_000_000)
            if pos.size:
                M = (G >= 0).astype(np.int8)
                X = np.where(G >= 0, G, 0).astype(np.int8)
                keep = L.filter_snps(X, M)
                pos, X, M = pos[keep], X[keep], M[keep]
            cnt, tot = L.binned_r2(pos, X, M, f) if pos.size else (np.zeros(L.NBIN, int), np.zeros(L.NBIN))
            for i in range(L.NBIN):
                fo.write(f"{label}\t{ch}\t{b}\t{f:.6f}\t{pos.size}\t{i}\t{cnt[i]}\t{tot[i]:.6f}\n")
            fo.flush()
            print(f"{label} {k}/{len(segs)} {ch}:{b} f={f:.2f} snps={pos.size} pairs={cnt.sum()}", flush=True)


if __name__ == "__main__":
    main()
