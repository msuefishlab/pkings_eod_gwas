#!/usr/bin/env python3
"""01_sim_null.py -- random-tract null from the local_ancestry simulator (code/11_Drift_vs_Selection/local_ancestry).

Regenerates simulated admixed fish (local_ancestry 05_simulate.simulate_hap: held-out donors, primary genome-wide
markers, r = 1.4 cM/Mb, donor redrawn only when ancestry changes) at each Apassa wobble's local_ancestry best-match
(T, m) cell plus the m = 0 cell; 50 fish per cell, seed 32032. Fits each simulated fish as local_ancestry
04_fit_hap.py fits real fish
(la_lib grid ML of (T, m) on both haplotypes; bootstrap B reduced to 50 because only the ML is used; posterior at
the ML) and saves per-marker diploid dosage and per-haplotype tract calls (posterior >= 0.5).

Writes <outdir>/sim_null.npz.
Usage (repo root): python3 code/14_Wobble_Analysis/apassa_ancestry/01_sim_null.py <la_dir> <outdir> [workers]
"""
import os, sys, csv
import numpy as np

LA_CODE = "code/11_Drift_vs_Selection/local_ancestry"
sys.path.insert(0, LA_CODE)
import la_lib as L
from importlib.machinery import SourceFileLoader
S = SourceFileLoader("sim05", f"{LA_CODE}/05_simulate.py").load_module()

la, out = sys.argv[1], sys.argv[2]
workers = int(sys.argv[3]) if len(sys.argv) > 3 else 12
os.makedirs(out, exist_ok=True)
NFISH, SEED = 50, 32032

cells = {(0, 0.0)}
for r in csv.DictReader(open(f"{la}/results/per_fish.tsv"), delimiter="\t"):
    if r["group"] in ("APA2019_Wobble", "APA2009_Wobble"):
        T = float(r["match_T"]); m = float(r["match_m"])
        cells.add((0, 0.0) if m == 0 else (int(T) if np.isfinite(T) else np.inf, m))
cells = sorted(cells, key=lambda x: (x[1], x[0]))
print("cells", cells, flush=True)

roles = S.read_roles(la)
M = L.load_markers(la, "primary")
gw = M["gw"].astype(bool); hs = M["hap_sample"]
colA = np.flatnonzero([("holdout_APA" in roles.get(s, ())) for s in hs])
colB = np.flatnonzero([("holdout_BAM" in roles.get(s, ())) for s in hs])
assert colA.size == 8 and colB.size == 20
chrom, pos = M["chrom"][gw], M["pos"][gw]
pBe, pTe = M["pBe"][gw], M["pTe"][gw]
HA, HB = M["H"][gw][:, colA], M["H"][gw][:, colB]
seg = L.segments(chrom, pos, L.excluded_regions())

rng = np.random.default_rng(SEED)
obs, meta = [], []
for T, m in cells:
    for rep in range(NFISH):
        for _ in range(2):
            al, _, _ = S.simulate_hap(rng, T, m, pos, seg, HA, HB)
            obs.append(al)
        meta.append((T, m, rep))
obs = np.stack(obs, 1)
print(f"simulated {len(meta)} fish", flush=True)
ones = np.ones(pos.size, bool)
Sh = L.grid_surfaces("hap", L.per_chrom_split(chrom, pos, obs, pBe, pTe, ones), workers)
Sh = Sh[..., 0::2] + Sh[..., 1::2]
fits = L.fit_from_surfaces(Sh, B=50); del Sh
lam = np.repeat([f["T"] * L.R_FIT for f in fits], 2); mm = np.repeat([f["m"] for f in fits], 2)
post = np.empty(obs.shape)
for c in range(1, 26):
    s = chrom == c
    post[s] = L.hap_posterior(obs[s], pos[s], pBe[s], pTe[s], lam, mm)
dose = (post[:, 0::2] + post[:, 1::2]).astype(np.float32)
callh = (post >= 0.5)
np.savez_compressed(f"{out}/sim_null.npz", dose=dose, call0=callh[:, 0::2], call1=callh[:, 1::2],
                    cell_T=np.array([x[0] for x in meta], float), cell_m=np.array([x[1] for x in meta]),
                    rep=np.array([x[2] for x in meta]), chrom=chrom, pos=pos,
                    T_hat=np.array([f["T"] for f in fits]), m_hat=np.array([f["m"] for f in fits]))
print("done", flush=True)
