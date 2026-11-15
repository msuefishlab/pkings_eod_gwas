#!/usr/bin/env python3
"""23b_sim_ld_decay.py -- simulated LD decay under one fitted A1w history at one r/mu (LD-decay test).

The demography is built exactly as in 04c_worker_msprime_agree.py (same params.json, same M3/M3f growth, same
migration naming). Each replicate is a 1-Mb callable-layout chunk with 67 BP and 47 TP diploids and
r = ratio x mu. The statistic is ld23_lib.binned_r2 on each population's dosages, with the same MAF filter as the
observed data (no missing data, f = 1).

Usage: 23b_sim_ld_decay.py <params.json> <model> <mu> <ratio> <n_reps> <out.tsv>
Output: one row per replicate x population x bin: model, mu, ratio, rep, pop, n_snps, bin, n_pairs, sum_r2.
"""
import json, os, sys
import numpy as np, msprime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ld23_lib as L

PARAMS, MODEL, MU, RATIO, NREP, OUT = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
P = json.load(open(PARAMS))[f"{MODEL}|{MU:g}"]
N_BP_IND, N_TP_IND = 67, 47
CHUNK = 1_000_000


def demography():
    # as 04c_worker_msprime_agree.py
    d = msprime.Demography()
    d.add_population(name="BP", initial_size=P["N_BP_end"])
    d.add_population(name="TP", initial_size=P["N_TP_end"])
    d.add_population(name="ANC", initial_size=P["N_anc"])
    if MODEL.rstrip("f") == "M3":
        T = P["T_gen"]
        d.populations[0].growth_rate = np.log(P["N_BP_end"] / P["N_BP_start"]) / T
        d.populations[1].growth_rate = np.log(P["N_TP_end"] / P["N_TP_start"]) / T
    if P["m_BP_from_TP_pergen"] > 0:
        d.set_migration_rate(source="BP", dest="TP", rate=P["m_BP_from_TP_pergen"])
    if P["m_TP_from_BP_pergen"] > 0:
        d.set_migration_rate(source="TP", dest="BP", rate=P["m_TP_from_BP_pergen"])
    d.add_population_split(time=P["T_gen"], derived=["BP", "TP"], ancestral="ANC")
    return d


DEM = demography()
seed0 = 20261030 + 1_000_000 * (sum(ord(c) for c in MODEL) % 97) + int(round(RATIO * 1000)) * 1000
with open(OUT, "w") as fo:
    fo.write("model\tmu\tratio\trep\tpop\tn_snps\tbin\tn_pairs\tsum_r2\n")
    for rep in range(1, NREP + 1):
        seed = seed0 + rep
        ts = msprime.sim_ancestry(samples={"BP": N_BP_IND, "TP": N_TP_IND}, demography=DEM,
                                  sequence_length=CHUNK, recombination_rate=RATIO * MU, ploidy=2,
                                  random_seed=seed)
        ts = msprime.sim_mutations(ts, rate=MU, model=msprime.BinaryMutationModel(),
                                   discrete_genome=False, random_seed=seed + 7)
        G = ts.genotype_matrix()                       # haploid, nodes in sample order (BP first)
        pos = ts.tables.sites.position
        for pop, a, b in (("BAM", 0, 2 * N_BP_IND), ("TP1", 2 * N_BP_IND, 2 * (N_BP_IND + N_TP_IND))):
            H = G[:, a:b]
            X = (H[:, 0::2] + H[:, 1::2]).astype(np.int8)  # diploid dosage: consecutive node pairs
            M = np.ones_like(X)
            keep = L.filter_snps(X, M)
            cnt, tot = L.binned_r2(pos[keep], X[keep], M[keep], 1.0)
            for i in range(L.NBIN):
                fo.write(f"{MODEL}\t{MU:g}\t{RATIO:g}\t{rep}\t{pop}\t{int(keep.sum())}\t{i}\t{cnt[i]}\t{tot[i]:.6f}\n")
        fo.flush()
