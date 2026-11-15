#!/usr/bin/env python3
"""
Simulate one chromosome for one demographic scenario with msprime.
Writes bgzip+tabix VCF to outdir/CHROM.vcf.gz

Recombination rate: r = 1e-8 per bp per gen (r/mu = 1).
NOTE: The smc++ estimate infers r/mu from data internally; we use r=1e-8
(not 6.22e-7) to keep the ARG tractable. With r/mu=62.2 the expected
tree count per chromosome exceeds ~10^9, making simulations infeasible.
The demographic recovery test is not biased by this choice.

Usage:
  python3 simulate_chrom.py --scenario null|geneflow|bottleneck
         --T_b INT  --chrom chrN  --length INT  --outdir DIR  --seed INT
"""
import sys, os, argparse, subprocess
import numpy as np
import msprime

MU     = 1e-8           # mutation rate per bp per gen
R      = 1e-8           # recombination rate per bp per gen (r/mu=1 for tractability)
N_ANCT = 47250          # ancient / plateau Ne (from pi-based estimate)
N_RECT = 2500           # recent Ne (bottleneck floor; midpoint of 1,500-3,000)
N_SAMP = 20             # diploid individuals
MASK   = 0.20           # fraction of sites to drop (approximates real masking)

def make_demography(scenario, T_b):
    if scenario == "null":
        d = msprime.Demography()
        d.add_population(name="POP", initial_size=N_ANCT)
        return d
    elif scenario == "geneflow":
        d = msprime.Demography()
        d.add_population(name="POP",   initial_size=N_ANCT)
        d.add_population(name="DONOR", initial_size=100_000)
        # 10 migrants/gen symmetric;  m = 10 / N_e
        m = 10.0 / N_ANCT
        d.set_symmetric_migration_rate(["POP", "DONOR"], m)
        return d
    else:  # bottleneck
        if T_b is None:
            raise ValueError("T_b required for bottleneck scenario")
        d = msprime.Demography()
        d.add_population(name="POP", initial_size=N_RECT)
        d.add_population_parameters_change(
            time=T_b, initial_size=N_ANCT, population="POP")
        return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", required=True,
                    choices=["null","geneflow","bottleneck"])
    ap.add_argument("--T_b",    type=int, default=None)
    ap.add_argument("--chrom",  required=True)
    ap.add_argument("--length", type=int, required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed",   type=int, default=42)
    args = ap.parse_args()

    vcf_gz  = os.path.join(args.outdir, f"{args.chrom}.vcf.gz")
    vcf_tbi = vcf_gz + ".tbi"

    if os.path.exists(vcf_tbi):
        print(f"[skip] {args.chrom}: already done", flush=True)
        sys.exit(0)

    os.makedirs(args.outdir, exist_ok=True)

    dem = make_demography(args.scenario, args.T_b)
    rng = np.random.default_rng(args.seed)

    # only sample from POP; DONOR is a ghost population for gene-flow scenario
    samples = {"POP": N_SAMP}

    print(f"[sim] {args.chrom}: len={args.length} scenario={args.scenario}"
          f" T_b={args.T_b} seed={args.seed}", flush=True)

    ts = msprime.sim_ancestry(
        samples=samples,
        demography=dem,
        sequence_length=args.length,
        recombination_rate=R,
        random_seed=args.seed,
    )
    ts = msprime.sim_mutations(ts, rate=MU, random_seed=args.seed + 1)

    # apply 20 % masking by dropping sites at random
    n_before = ts.num_sites
    if n_before > 0:
        keep         = rng.random(n_before) > MASK
        to_delete    = np.where(~keep)[0].tolist()
        ts           = ts.delete_sites(to_delete)

    print(f"[sim] {args.chrom}: {ts.num_sites} sites kept (from {n_before})",
          flush=True)

    # write VCF
    vcf_tmp = vcf_gz.replace(".vcf.gz", ".vcf")
    with open(vcf_tmp, "w") as fh:
        ts.write_vcf(fh, contig_id=args.chrom,
                     position_transform=lambda x: 1 + x)

    # bgzip + tabix (requires HTSlib module loaded in SLURM environment)
    subprocess.run(["bgzip", "-f", vcf_tmp], check=True)
    subprocess.run(["tabix", "-p", "vcf", vcf_gz], check=True)

    print(f"[done] {args.chrom} -> {vcf_gz}", flush=True)

if __name__ == "__main__":
    main()
