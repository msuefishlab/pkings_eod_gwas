#!/usr/bin/env python3
"""
build_masks.py — per-population callability mask for smc++

This is only the call-rate component of the original masks. It cannot see positions
that are absent from the all-sites VCF (indels, multiallelic and spanning-deletion
records and non-PASS invariant blocks are dropped by 03_VCF_QC_Filter; 14.15% of
chr1-25), and vcf2smc reads an unmasked absent position as homozygous reference.
build_masks_v2.sb therefore adds every VCF-absent interval to the original masks;
use only those masks_v2 with smc++. generate_awk.py + compute_callrate.awk implement
the same call-rate rule as a per-chromosome awk scan.

Single-pass through bcftools query output, computing per-site call rate
for each of 8 populations simultaneously.

Usage:
    bcftools query -f '%CHROM\t%POS[\t%GT]\n' \
        -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25 \
        <vcf.gz> |
    python3 build_masks.py \
        --vcf <vcf.gz> \
        --samples-dir <dir> \
        --outdir <dir> \
        --pops BAM BIK BAVA BIR APA_2009 BEN DOV DOG

Outputs: <outdir>/<POP>.lowcall.bed (0-based half-open BED;
         sites where per-population call rate < min_callrate)
"""

import sys
import os
import subprocess
import argparse


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--vcf', required=True,
                   help='Path to VCF (for sample header extraction with bcftools query -l)')
    p.add_argument('--samples-dir', required=True,
                   help='Directory containing <POP>.txt sample lists')
    p.add_argument('--outdir', required=True,
                   help='Output directory for <POP>.lowcall.bed files')
    p.add_argument('--pops', nargs='+', required=True,
                   help='Population names (must match <POP>.txt files in samples-dir)')
    p.add_argument('--min-callrate', type=float, default=0.8,
                   help='Minimum acceptable call rate (default: 0.8)')
    return p.parse_args()


def load_pop_samples(samples_dir, pops):
    """Return dict: pop -> frozenset of sample IDs."""
    pop_samples = {}
    for pop in pops:
        path = os.path.join(samples_dir, f'{pop}.txt')
        with open(path) as fh:
            samples = [
                line.strip() for line in fh
                if line.strip() and not line.startswith('#')
            ]
        pop_samples[pop] = frozenset(samples)
        print(f"[build_masks] {pop}: {len(samples)} samples loaded from {path}",
              file=sys.stderr, flush=True)
    return pop_samples


def get_vcf_sample_order(vcf_path):
    """Return list of sample IDs in VCF column order using bcftools query -l."""
    result = subprocess.run(
        ['bcftools', 'query', '-l', vcf_path],
        capture_output=True, text=True, check=True
    )
    samples = result.stdout.strip().split('\n')
    print(f"[build_masks] VCF has {len(samples)} samples", file=sys.stderr, flush=True)
    return samples


def build_pop_indices(vcf_samples, pop_samples):
    """Return dict: pop -> list of 0-based column indices in vcf_samples."""
    vcf_sample_set = {s: i for i, s in enumerate(vcf_samples)}
    pop_indices = {}
    for pop, samples in pop_samples.items():
        indices = []
        missing = []
        for s in samples:
            if s in vcf_sample_set:
                indices.append(vcf_sample_set[s])
            else:
                missing.append(s)
        if missing:
            print(f"[build_masks] WARNING {pop}: {len(missing)} samples not in VCF: "
                  f"{missing[:5]}{'...' if len(missing) > 5 else ''}",
                  file=sys.stderr, flush=True)
        pop_indices[pop] = indices
        print(f"[build_masks] {pop}: {len(indices)} samples mapped to VCF columns",
              file=sys.stderr, flush=True)
    return pop_indices


def is_missing(gt_str):
    """Return True if genotype call is missing (./. or . or .|. with . allele)."""
    # Normalise phased/unphased
    alleles = gt_str.replace('|', '/').split('/')
    return alleles[0] == '.'


def main():
    args = parse_args()

    pop_samples = load_pop_samples(args.samples_dir, args.pops)
    vcf_samples = get_vcf_sample_order(args.vcf)
    pop_indices = build_pop_indices(vcf_samples, pop_samples)

    os.makedirs(args.outdir, exist_ok=True)

    # Open one output file per population
    outfiles = {
        pop: open(os.path.join(args.outdir, f'{pop}.lowcall.bed'), 'w')
        for pop in args.pops
    }

    target_chroms = {f'chr{i}' for i in range(1, 26)}

    n_lines = 0
    n_written = {pop: 0 for pop in args.pops}

    for raw_line in sys.stdin:
        parts = raw_line.rstrip('\n').split('\t')
        chrom = parts[0]
        if chrom not in target_chroms:
            continue
        pos = int(parts[1])   # 1-based VCF position
        gts = parts[2:]       # one GT string per VCF sample, in VCF column order

        for pop, indices in pop_indices.items():
            n_total = len(indices)
            if n_total == 0:
                continue
            n_missing = sum(1 for i in indices if is_missing(gts[i]))
            call_rate = (n_total - n_missing) / n_total
            if call_rate < args.min_callrate:
                # 0-based half-open BED: [pos-1, pos)
                outfiles[pop].write(f"{chrom}\t{pos - 1}\t{pos}\n")
                n_written[pop] += 1

        n_lines += 1
        if n_lines % 500_000 == 0:
            report = ', '.join(f'{p}:{n_written[p]:,}' for p in args.pops)
            print(f"[build_masks] {n_lines:,} sites processed — low-call: {report}",
                  file=sys.stderr, flush=True)

    for fh in outfiles.values():
        fh.close()

    print(f"[build_masks] Done. {n_lines:,} sites processed.", file=sys.stderr, flush=True)
    for pop in args.pops:
        print(f"[build_masks]   {pop}: {n_written[pop]:,} low-callrate sites written",
              file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
