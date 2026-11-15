#!/usr/bin/env python3
"""
17_xpnsl.py
XP-nSL (cross-population nSL) for one BP-TP contrast, one chromosome at a time (the XP-nSL
sensitivity row of Table S5).

XP-nSL compares haplotype homozygosity between two populations at each site, measured in segregating
sites rather than base pairs, so it depends less on the local recombination rate than Rsb or XP-EHH.
It is used as an alternative to Rsb for the haplotype class of the sweep-evidence rule.

Input is the same as the Rsb scan (scan_ehh_integrated.R): the phased pair VCF, restricted to sites
whose INFO/AA matches REF or ALT (rehh's polarize_vcf = TRUE drops the rest), filtered to pooled
MAF >= 0.05, then split into BP and TP with the scan's own sample lists. Sites with any missing
haplotype are dropped. XP-nSL is not allele-specific, so the polarization itself does not change a
score; the restriction only keeps the site set identical to Rsb's.

Raw (unstandardized) scores are written; 17b_xpnsl_merge.py standardizes them genome-wide. Positive
means longer haplotypes in BP.

Usage (inside dadi_tools.sif):
  python3 17_xpnsl.py --vcf <pair VCF> --bp-list <BP.txt> --tp-list <TP.txt> --chrom chr6 --out <out.tsv.gz>
"""

import argparse
import gzip
import sys

import numpy as np

try:
    import allel
except ImportError:
    sys.exit("scikit-allel is required (run this inside dadi_tools.sif)")


def read_ids(fn):
    with open(fn) as fh:
        return [line.split()[0] for line in fh if line.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="phased pair VCF for one chromosome")
    ap.add_argument("--bp-list", required=True)
    ap.add_argument("--tp-list", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--out", required=True, help="output .tsv.gz (chrom, pos, raw)")
    ap.add_argument("--min-maf", type=float, default=0.05)
    args = ap.parse_args()

    cb = allel.read_vcf(args.vcf, fields=["samples", "variants/POS", "variants/REF", "variants/ALT",
                                          "variants/AA", "calldata/GT"])
    if cb is None or "calldata/GT" not in cb:
        sys.exit(f"no genotypes read from {args.vcf}")
    samples = list(cb["samples"])
    idx = {s: i for i, s in enumerate(samples)}
    bp = [idx[s] for s in read_ids(args.bp_list) if s in idx]
    tp = [idx[s] for s in read_ids(args.tp_list) if s in idx]
    print(f"{args.chrom}: {len(samples)} samples in VCF; BP {len(bp)}, TP {len(tp)} matched")
    if len(bp) < 5 or len(tp) < 5:
        sys.exit("too few samples matched the lists")

    pos = cb["variants/POS"]
    gt = allel.GenotypeArray(cb["calldata/GT"])
    h_all = gt.take(bp + tp, axis=1).to_haplotypes()
    hv = np.asarray(h_all)
    aa = np.array([str(a) if a is not None else "." for a in cb["variants/AA"]])
    polarized = (aa == cb["variants/REF"]) | (aa == cb["variants/ALT"][:, 0])
    complete = (hv >= 0).all(axis=1)
    freq = np.where(complete, hv.clip(0).sum(axis=1) / hv.shape[1], 0.0)
    keep = polarized & complete & (freq >= args.min_maf) & (freq <= 1 - args.min_maf)
    print(f"  {len(pos)} variants; {polarized.sum()} polarized; {complete.sum()} complete; "
          f"{keep.sum()} kept (polarized, complete, pooled MAF >= {args.min_maf})")
    if keep.sum() < 100:
        sys.exit("too few sites to compute XP-nSL")

    hk = hv[keep]
    n_bp = 2 * len(bp)
    h1 = allel.HaplotypeArray(hk[:, :n_bp])
    h2 = allel.HaplotypeArray(hk[:, n_bp:])
    score = allel.xpnsl(h1, h2, use_threads=True)

    ok = np.isfinite(score)
    print(f"  {ok.sum()} finite raw scores")
    with gzip.open(args.out, "wt") as out:
        out.write("chrom\tpos\traw\n")
        for p, v in zip(pos[keep][ok], score[ok]):
            out.write(f"{args.chrom}\t{p}\t{v:.6f}\n")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
