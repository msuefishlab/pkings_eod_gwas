#!/usr/bin/env python3
"""One task: simulate 5-Mb chunks of the fitted A1w history until a wall-clock budget is spent.
Writes one row per chunk (marginal counts) and a Counter of (d_BP,d_TP,het_BP,het_TP) for delta>=0.70.
Usage: 04c_worker_msprime_agree.py <model> <mu> <task_id> <budget_s> <outdir>
"""
import sys, time, json, collections, os
import numpy as np, msprime

MODEL, MU, TASK, BUDGET, OUTDIR = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
PAR = json.load(open(os.path.join(OUTDIR, "params.json")))
P = PAR[f"{MODEL}|{MU:g}"]
R_RATE = PAR["r_rate"]; CHUNK = int(PAR["chunk_len"])
N_BP_IND, N_TP_IND = 67, 47
DELTA_MIN = 0.70
THRESH = {"peak4": PAR["delta_peak4"], "peak5": PAR["delta_peak5"]}
EPS = 1e-9

def demography():
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
rows = []
posrows = []
joint = collections.Counter()
t0 = time.time()
chunk = 0
seed0 = 1_000_000 * (1 + sum(ord(ch) for ch in MODEL) % 97) + 100_000 * (1 if MU < 7e-9 else 2) + 1000 * TASK
while time.time() - t0 < BUDGET:
    chunk += 1
    seed = seed0 + chunk
    ts = msprime.sim_ancestry(samples={"BP": N_BP_IND, "TP": N_TP_IND}, demography=DEM,
                             sequence_length=CHUNK, recombination_rate=R_RATE, ploidy=2,
                             random_seed=seed)
    ts = msprime.sim_mutations(ts, rate=MU, model=msprime.BinaryMutationModel(),
                               discrete_genome=False, random_seed=seed + 7)
    if chunk == 1:
        samp = ts.samples()
        popn = [ts.node(j).population for j in samp]
        bpid = ts.population(0).id
        assert len(samp) == 228 and all(p == bpid for p in popn[:134]) and all(p != bpid for p in popn[134:]), "sample order"
        ind_nodes = [tuple(ts.individual(i).nodes) for i in range(5)]
        assert all(n[1] == n[0] + 1 for n in ind_nodes), f"node pairing {ind_nodes}"
    G = ts.genotype_matrix()                      # (n_sites, 228) haploid, nodes in sample order
    if G.shape[0] == 0:
        continue
    bp = G[:, :2 * N_BP_IND]; tp = G[:, 2 * N_BP_IND:]
    d_bp = bp.sum(1).astype(np.int32); d_tp = tp.sum(1).astype(np.int32)
    # individuals are consecutive node pairs (ploidy=2)
    het_bp = (bp[:, 0::2] != bp[:, 1::2]).sum(1).astype(np.int32)
    het_tp = (tp[:, 0::2] != tp[:, 1::2]).sum(1).astype(np.int32)
    f_bp = d_bp / (2 * N_BP_IND); f_tp = d_tp / (2 * N_TP_IND)
    delta = f_bp - f_tp
    row = dict(model=MODEL, mu=MU, task=TASK, chunk=chunk, seed=seed, L=CHUNK, n_sites=int(G.shape[0]))
    for nm, thr in THRESH.items():
        row[f"n_{nm}_der"] = int((delta >= thr - EPS).sum())          # mutant allele BP-biased
        row[f"n_{nm}_anc"] = int((-delta >= thr - EPS).sum())         # ancestral allele BP-biased
        row[f"n_{nm}_either"] = row[f"n_{nm}_der"] + row[f"n_{nm}_anc"]
    rows.append(row)
    keep = np.flatnonzero(np.abs(delta) >= DELTA_MIN - EPS)
    for i in keep:
        joint[(int(d_bp[i]), int(d_tp[i]), int(het_bp[i]), int(het_tp[i]))] += 1
    thr5 = THRESH["peak5"]
    hit = np.flatnonzero(delta >= thr5 - EPS)
    if hit.size:
        pos = ts.tables.sites.position[hit]
        for i, p_ in zip(hit, pos):
            posrows.append((chunk, float(p_), 1 if delta[i] >= THRESH["peak4"] - EPS else 0))

tag = f"{MODEL}_{MU:g}_{TASK}"
import csv
with open(os.path.join(OUTDIR, "chunks", f"chunks_{tag}.tsv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
with open(os.path.join(OUTDIR, "joint", f"joint_{tag}.tsv"), "w") as fh:
    fh.write("model\tmu\ttask\td_BP\td_TP\thet_BP\thet_TP\tcount\n")
    for (a, b, cc, dd), n in sorted(joint.items()):
        fh.write(f"{MODEL}\t{MU:g}\t{TASK}\t{a}\t{b}\t{cc}\t{dd}\t{n}\n")
with open(os.path.join(OUTDIR, "pos", f"pos_{tag}.tsv"), "w") as fh:
    fh.write("model\tmu\ttask\tchunk\tpos\tge_peak4\n")
    for cnum, p_, g4 in posrows:
        fh.write(f"{MODEL}\t{MU:g}\t{TASK}\t{cnum}\t{p_:.0f}\t{g4}\n")
print(f"{tag}: {chunk-1} chunks, {sum(r['n_sites'] for r in rows)} sites, {time.time()-t0:.0f}s")
