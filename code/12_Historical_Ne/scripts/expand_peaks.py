#!/usr/bin/env python3
"""
expand_peaks.py — expand GWAS peaks BED ± N bp and add chr prefix.

Usage:
    python3 expand_peaks.py <peaks.bed> <genome.fai> [--expand 1000000]

Outputs a BED file (stdout or --out) with:
  - chr prefix added to column 1 if absent
  - each peak expanded by --expand bp on each side
  - clipped to chromosome bounds from the FAI
"""

import sys
import os
import argparse


def main():
    p = argparse.ArgumentParser()
    p.add_argument('peaks_bed', help='Input peaks BED (tab-separated)')
    p.add_argument('fai', help='Reference .fai file for chromosome lengths')
    p.add_argument('--expand', type=int, default=1_000_000,
                   help='Expand each peak by this many bp on each side (default 1000000)')
    p.add_argument('--out', default='-', help='Output BED file (default: stdout)')
    args = p.parse_args()

    # Read chromosome lengths from FAI
    chrom_lengths = {}
    with open(args.fai) as fh:
        for line in fh:
            parts = line.strip().split('\t')
            chrom_lengths[parts[0]] = int(parts[1])

    out = open(args.out, 'w') if args.out != '-' else sys.stdout

    n = 0
    with open(args.peaks_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom = parts[0]
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom
            start = max(0, int(parts[1]) - args.expand)
            end = int(parts[2]) + args.expand
            if chrom in chrom_lengths:
                end = min(end, chrom_lengths[chrom])
            out.write(f"{chrom}\t{start}\t{end}\n")
            n += 1

    if args.out != '-':
        out.close()

    print(f"[expand_peaks] {n} peaks written (each expanded by ±{args.expand:,} bp)",
          file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
