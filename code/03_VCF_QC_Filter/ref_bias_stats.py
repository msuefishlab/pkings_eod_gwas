#!/usr/bin/env python3
"""Per-sample reference-bias statistics for one scaffold of the raw all-sites VCF.

For each region class in the BED from ref_bias_make_regions.py, and each sample:
  depth pass (invariant + SNP sites): n_sites, sum_dp, n_miss (DP missing or 0)
  allele-balance pass (biallelic SNPs, het calls, DP >= 10): n_het, sum_frac (ref reads / (ref+alt)),
  sum_ref, sum_alt
A BP reference that TP reads map to poorly would show as lower relative depth, more missing calls and
ref-skewed heterozygotes in TP fish, most of all in the peaks and the highly variable windows.

Usage: python3 ref_bias_stats.py <vcf.gz> <regions.bed> <scaffold> <out.tsv>
"""
import sys, subprocess, tempfile, os
import numpy as np, pandas as pd

vcf, bed, scaf, out = sys.argv[1:5]
CHUNK = 100000

reg = pd.read_csv(bed, sep="\t", header=None, names=["scaf", "start", "end", "label", "chr"])
reg = reg[reg.scaf == scaf].sort_values("start").reset_index(drop=True)
starts, ends = reg.start.to_numpy(), reg.end.to_numpy()
labels = sorted(reg.label.unique())
lab_idx = reg.label.map({l: i for i, l in enumerate(labels)}).to_numpy()

tmp = tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False)
reg[["scaf", "start", "end"]].to_csv(tmp.name, sep="\t", header=False, index=False)
tmp.close()

samples = subprocess.run(["bcftools", "query", "-l", vcf], capture_output=True, text=True, check=True).stdout.split()
ns, nl = len(samples), len(labels)

def label_of(pos):
    """Region-class index for 1-based positions (regions are non-overlapping)."""
    i = np.searchsorted(starts, pos - 1, side="right") - 1
    ok = (i >= 0) & (pos - 1 < ends[np.clip(i, 0, None)])
    return np.where(ok, lab_idx[np.clip(i, 0, None)], -1)

def stream(args):
    return subprocess.Popen(["bcftools", "query", "-R", tmp.name] + args + [vcf], stdout=subprocess.PIPE)

# depth pass
n_sites = np.zeros(nl)
sum_dp = np.zeros((nl, ns))
n_miss = np.zeros((nl, ns))
p = stream(["-i", 'TYPE="ref" || TYPE="snp"', "-f", "%POS[\t%DP]\n"])
for ch in pd.read_csv(p.stdout, sep="\t", header=None, na_values=".", chunksize=CHUNK):
    pos = ch.iloc[:, 0].to_numpy()
    dp = ch.iloc[:, 1:].to_numpy(dtype=float)
    li = label_of(pos)
    for k in range(nl):
        m = li == k
        if m.any():
            d = dp[m]
            n_sites[k] += m.sum()
            sum_dp[k] += np.nansum(d, axis=0)
            n_miss[k] += (np.isnan(d) | (d == 0)).sum(axis=0)
assert p.wait() == 0, "depth pass failed"

# allele-balance pass. bcftools query -i applies FORMAT tests per site, not per sample, so het calls
# are picked out here.
n_het = np.zeros((nl, ns))
sum_frac = np.zeros((nl, ns))
sum_ref = np.zeros((nl, ns))
sum_alt = np.zeros((nl, ns))
HET = ["0/1", "0|1", "1|0"]
p = stream(["-i", 'TYPE="snp" && N_ALT=1', "-f", "%POS[\t%GT\t%DP\t%AD{0}\t%AD{1}]\n"])
for ch in pd.read_csv(p.stdout, sep="\t", header=None, na_values=".", chunksize=CHUNK):
    li = label_of(ch.iloc[:, 0].to_numpy())
    v = ch.iloc[:, 1:]
    gt = v.iloc[:, 0::4].to_numpy(dtype=str)
    dp = v.iloc[:, 1::4].to_numpy(dtype=float)
    r = v.iloc[:, 2::4].to_numpy(dtype=float)
    a = v.iloc[:, 3::4].to_numpy(dtype=float)
    ok = np.isin(gt, HET) & (dp >= 10) & ((r + a) > 0) & (li >= 0)[:, None]
    frac = np.where(ok, r / np.where(ok, r + a, 1), 0)
    for k in range(nl):
        m = li == k
        if m.any():
            n_het[k] += ok[m].sum(axis=0)
            sum_frac[k] += frac[m].sum(axis=0)
            sum_ref[k] += np.where(ok[m], r[m], 0).sum(axis=0)
            sum_alt[k] += np.where(ok[m], a[m], 0).sum(axis=0)
assert p.wait() == 0, "allele-balance pass failed"
os.unlink(tmp.name)

with open(out, "w") as o:
    o.write("scaf\tlabel\tsample\tn_sites\tsum_dp\tn_miss\tn_het\tsum_frac\tsum_ref\tsum_alt\n")
    for k, l in enumerate(labels):
        for j, s in enumerate(samples):
            o.write(f"{scaf}\t{l}\t{s}\t{int(n_sites[k])}\t{int(sum_dp[k, j])}\t{int(n_miss[k, j])}\t"
                    f"{int(n_het[k, j])}\t{sum_frac[k, j]:.4f}\t{int(sum_ref[k, j])}\t{int(sum_alt[k, j])}\n")
