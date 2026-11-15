#!/usr/bin/env python3
"""lrt_ext.py -- Godambe-adjusted LRTs for information only: CB vs M0, CBT vs CB, CBTM vs CBT, in the
unfolded and the folded analysis. Mirrors the LRT block of the drift test's 08_collect_fits.py: D_raw = 2 (LL_complex - LL_simple);
dadi.Godambe.LRT_adjust(func_ex of the complex model, pts_l, the 100 x 1-Mb bootstrap spectra, p0, data,
nested_indices, multinom=True); D_adj = adj * D_raw; p from the chi-square mixture for parameters on a boundary
(Self & Liang 1987), as 08: one parameter 0.5 chi2_1, two parameters 0.5 chi2_1 + 0.25 chi2_2.
p0 is the complex model's parameter vector at the simple model's optimum:
  CB vs M0:     (nu_BP, nu_TP from M0; T_pre = T_M0 - T_c*, T_c*; r_BP = 1; p_misid from M0), T_c* = CB's fitted T_c
                if it is below T_M0 - 1e-5, else T_M0 / 2 (T_c is not identified under the null)
  CBT vs CB:    CB's fit, r_TP = 1
  CBTM vs CBT:  CBT's fit, m = 0, 0, f_mig = CBTM's fitted f_mig (not identified under the null)
Nested parameters: r_BP = 1 and r_TP = 1 sit on the upper bound; m = 0 on the lower bound. Because T_c (CB vs M0) and
f_mig (CBTM vs CBT) are not identified under the null, the mixture p is nominal (Davies problem); information only.
Usage: FD_OUT=<fits root> python3 lrt_ext.py --pairs A1w --out lrt.tsv
"""
import argparse
import csv
import importlib.util
import multiprocessing as mp
import os
import pickle
import sys

import dadi
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ext_models  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--pairs", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
OUT = os.path.join(os.getcwd(), os.environ["FD_OUT"])
spec = importlib.util.spec_from_file_location("fitmod", os.path.join(HERE, "..", "..", "fixed_differences", "07_fit_models.py"))
fitmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fitmod)
ext_models.register(fitmod)
best = {(r["pair"], r["model"]): r for r in csv.DictReader(open(os.path.join(OUT, "fits/best_fits.tsv")), delimiter="\t")}


def P(r):
    return {k: float(x) for k, x in (t.split("=") for t in r["params"].split(";"))}


def one(job):
    pair, simple, comp = job
    if (pair, simple) not in best or (pair, comp) not in best:
        return None
    f = "f" if comp.endswith("f") else ""
    ps, pc = P(best[(pair, simple)]), P(best[(pair, comp)])
    b = comp.rstrip("f")
    note = ""
    if b == "CB":
        T = ps["T"]
        tc = pc["T_c"] if pc["T_c"] < T - 1e-5 else T / 2
        note = f"T_c*={tc:.6g}" + ("" if tc == pc["T_c"] else " (T/2)")
        q = {"nu_BP": ps["nu_BP"], "nu_TP": ps["nu_TP"], "T_pre": T - tc, "T_c": tc, "r_BP": 1.0}
        nested = [4]
    elif b == "CBT":
        q = {k: ps[k] for k in ext_models.PNAMES["CB"]}
        q["r_TP"] = 1.0
        nested = [5]
    else:
        q = {k: ps[k] for k in ext_models.PNAMES["CBT"]}
        q.update({"m_BP_from_TP": 0.0, "m_TP_from_BP": 0.0, "f_mig": pc["f_mig"]})
        note = f"f_mig={pc['f_mig']:.6g}"
        nested = [6, 7]
    p0 = [q[k] for k in ext_models.PNAMES[b]] + ([] if f else [ps["p_misid"]])
    pk = pickle.load(open(os.path.join(OUT, f"sfs/{pair}.pkl"), "rb"))
    data = fitmod.data_of(pk, comp)
    func_ex, pts_l, pnames, _, _, _ = fitmod.model_setup(comp, data.sample_sizes)
    boots = fitmod.bootstrap_spectra(pk, block_mb=1, folded=fitmod.is_folded(comp))
    d_raw = 2 * (float(best[(pair, comp)]["LL"]) - float(best[(pair, simple)]["LL"]))
    try:
        adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, boots, p0, data, nested_indices=nested, multinom=True)
        d_adj = adj * d_raw
        p = 0.5 * chi2.sf(d_adj, 1) if len(nested) == 1 else 0.5 * chi2.sf(d_adj, 1) + 0.25 * chi2.sf(d_adj, 2)
        err = ""
    except Exception as e:
        adj = d_adj = p = float("nan")
        err = str(e)
    print(pair, simple, comp, f"D_raw {d_raw:.3f} adj {adj:.5g} D_adj {d_adj:.3f} p {p:.4g}", flush=True)
    return dict(pair=pair, analysis="folded" if f else "unfolded", simple=simple, complex=comp,
                nested=",".join(pnames[i] for i in nested), LL_simple=best[(pair, simple)]["LL"],
                LL_complex=best[(pair, comp)]["LL"], D_raw=f"{d_raw:.3f}", adj_factor=f"{adj:.6g}",
                D_adj=f"{d_adj:.3f}", p_adj_mixture=f"{p:.4g}", p0=",".join(f"{v:.6g}" for v in p0), note=note,
                error=err)


jobs = [(p, s + f, c + f) for p in a.pairs.split(",") for f in ("", "f")
        for s, c in (("M0", "CB"), ("CB", "CBT"), ("CBT", "CBTM"))]
n = min(len(jobs), int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
with mp.get_context("fork").Pool(max(n, 1)) as pool:
    rows = [r for r in pool.map(one, jobs, chunksize=1) if r]
cols = list(rows[0]) if rows else ["pair"]
with open(a.out + ".tmp", "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in rows:
        o.write("\t".join(str(r.get(k, "")) for k in cols) + "\n")
os.replace(a.out + ".tmp", a.out)
