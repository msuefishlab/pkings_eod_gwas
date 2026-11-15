#!/usr/bin/env python3
"""Collate smc++ results (masks_v2 pipeline). Usage: collate.py <work_tree> <out_dir>.
t = knot*2*N0, Ne = exp(y)*N0, piecewise-constant interpolation on 200 log-spaced
points 10..2e6 gen, 95% bootstrap percentile bands; split T = t_split*2*N0 of the
pop1 marginal; r/mu reported raw (rho/theta) and divided by the window size w = 100
(column r_per_mu_corrected_w100 of rho_mu_table.tsv)."""
import os, sys, json
import numpy as np
V2 = sys.argv[1]; OUT = sys.argv[2]; os.makedirs(OUT, exist_ok=True)
POPS = ["BAM", "BIK", "BAVA", "BIR", "APA_2009", "BEN", "DOV", "DOG"]
LIN = {"BAM": "BP", "BIK": "BP", "BAVA": "BP", "BIR": "BP", "APA_2009": "TP", "BEN": "TP", "DOV": "TP", "DOG": "TP"}
MU = [("1e-8", "1e8"), ("4e-9", "4e9")]
PAIRS = [("BAM", "APA_2009"), ("BIK", "DOV"), ("BAVA", "DOV"), ("BIR", "DOV"), ("APA_2009", "BEN"), ("DOV", "DOG")]
W_CORR = 100.0
t_ref = np.logspace(1, np.log10(2e6), 200)

def load(p):
    with open(p) as f: d = json.load(f)
    m = d["model"]; N0 = float(m["N0"])
    return np.array(m["knots"], float) * 2 * N0, np.exp(np.array(m["y"], float)) * N0, d.get("rho"), d.get("theta"), N0, d
def interp(t, tk, nk):
    i = np.clip(np.searchsorted(tk, t, side="right") - 1, 0, len(nk) - 1); return nk[i]

rows, rho_rows, status = [], [], []
for pop in POPS:
    for mu, lab in MU:
        pj = f"{V2}/estimates/{pop}_mu{lab}/model.final.json"
        if not os.path.exists(pj): status.append(f"MISSING point {pop} {lab}"); continue
        tk, nk, rho, theta, N0, _ = load(pj); Np = interp(t_ref, tk, nk)
        boots = []
        for K in range(1, 11):
            bj = f"{V2}/bootstraps/estimates/boot_{K}/{pop}_mu{lab}/model.final.json"
            if os.path.exists(bj):
                bt, bn, *_ = load(bj); boots.append(interp(t_ref, bt, bn))
            else: status.append(f"MISSING boot {K} {pop} {lab}")
        lo = np.percentile(boots, 2.5, axis=0) if boots else np.full(200, np.nan)
        hi = np.percentile(boots, 97.5, axis=0) if boots else np.full(200, np.nan)
        for i, t in enumerate(t_ref):
            rows.append((pop, lab, t, Np[i], lo[i], hi[i], len(boots)))
        os.makedirs(f"{OUT}/Ne_tables", exist_ok=True)   # same format as the original per-pop tables
        with open(f"{OUT}/Ne_tables/{pop}_mu{lab}_Ne_table.tsv", "w") as f:
            f.write("t_gen\tNe\tNe_lo_95\tNe_hi_95\tt_yr_1y\tt_yr_2y\n")
            for i, t in enumerate(t_ref):
                f.write(f"{t:.2f}\t{Np[i]:.2f}\t{lo[i]:.2f}\t{hi[i]:.2f}\t{t*1:.2f}\t{t*2:.2f}\n")
        rmu = float(rho) / float(theta) if rho is not None and theta else float("nan")
        rho_rows.append((LIN[pop], pop, lab, rho, theta, rmu, rmu / W_CORR, N0, len(boots)))

with open(f"{OUT}/Ne_all_pops_both_mu.tsv", "w") as f:
    f.write("pop\tmu\tt_gen\tNe\tNe_lo_95\tNe_hi_95\tt_yr_1y\tt_yr_2y\tn_boot\n")
    for p, l, t, n, a, b, k in rows:
        f.write(f"{p}\t{l}\t{t:.2f}\t{n:.2f}\t{a:.2f}\t{b:.2f}\t{t:.2f}\t{2*t:.2f}\t{k}\n")
with open(f"{OUT}/rho_mu_table.tsv", "w") as f:
    f.write("lineage\tpop\tmu\trho\ttheta\tr_per_mu_reported\tr_per_mu_corrected_w100\tN0\tn_boot\n")
    for r in rho_rows: f.write("\t".join(str(x) for x in r) + "\n")

# Both split runs go in one table. mode = flags used (folded adds --polarization-error 0.5,
# which smc++ 1.15.4 split ignores), so the two rows per pair x mu are replicate fits.
with open(f"{OUT}/split_table.tsv", "w") as f:
    f.write("pair\tpop1\tpop2\tmu\tmode\tt_split_scaled\tN0_pop1\tT_split_gen\tT_split_yr_1y\tT_split_yr_2y\n")
    for sub, tag in [("split", "unfolded"), ("split_folded", "folded")]:
        for p1, p2 in PAIRS:
            for _, lab in MU:
                pair = f"{p1}-{p2}"; sj = f"{V2}/{sub}/{pair}_mu{lab}/model.final.json"
                mj = f"{V2}/estimates/{p1}_mu{lab}/model.final.json"
                if not (os.path.exists(sj) and os.path.exists(mj)):
                    status.append(f"MISSING {sub} {pair} {lab}"); continue
                sd = json.load(open(sj)); N0 = json.load(open(mj))["model"]["N0"]
                ts = None
                for d in (sd, sd.get("model", {})):
                    for k in ("split", "t_split", "T_split", "split_t", "split_time"):
                        if ts is None and k in d: ts = float(d[k])
                if ts is None: status.append(f"NO split field {sj}"); continue
                T = ts * 2 * N0
                f.write(f"{pair}\t{p1}\t{p2}\tmu{lab}\t{tag}\t{ts}\t{N0}\t{round(T)}\t{round(T)}\t{round(2*T)}\n")

with open(f"{OUT}/collation_status.txt", "w") as f:
    f.write("\n".join(status) + ("\n" if status else "ALL PRESENT\n"))
print(f"rows={len(rows)} rho_rows={len(rho_rows)} problems={len(status)}")
