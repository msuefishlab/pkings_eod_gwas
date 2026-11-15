#!/usr/bin/env python3
"""gate_ext.py -- Stage 2b fit status per pair x model (no index-SNP / Peak statistic is read or computed).
Gate: the per-model admissibility block of the drift test's 09_tail.py (theta-scaled expected spectrum,
Delta ladder 0.5..1.0 over the fitting mask, gate = most extreme threshold with E >= 10, adequate iff 1/3 <= O/E <= 3,
admissible = converged and adequate and (folded or p_misid <= 0.05)). It does not open tail/index_units.tsv.
Stage 2b additions, for CB / CBT / CBTM (+ f) only:
  nested_ok_ext  nested-start rule of the drift test along ext_models.SUBMODELS (same folding; wide run: also the same model at primary
                 bounds, --refdir): best LL below any sub-model's best LL -> not converged
  converged      = 08's Godambe-adjusted convergence AND nested_ok_ext
  bounds         box of the best start (fits/bounds/*.json sidecar); a parameter within 1% (ratio) of a positive bound,
                 or < 1e-3 at a zero lower bound, is "at bound"
  T_c window     exact x = theta_opt*T_c/(2L) against the generation window (tolerance 1% in x)
  SEs            GIM recomputed with the call 08 makes (same bootstrap spectra, same parameters); cov = GIM^-1 over
                 (free params, theta); SEs of derived (nu_BP2, nu_TP2, T_s, T_m) and physical quantities (both mu arms)
                 by the delta method (theta uncertainty included); the direct SEs are checked against 08's
Usage: FD_OUT=<fits root> python3 gate_ext.py --pairs A1w --models M0,...,CBTMf --bounds primary --out f.tsv
"""
import argparse
import csv
import glob
import importlib.util
import json
import multiprocessing as mp
import os
import pickle
import sys

import numpy as np
import dadi

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ext_models  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--pairs", required=True)
ap.add_argument("--models", required=True)
ap.add_argument("--bounds", default="primary")
ap.add_argument("--refdir", default=None, help="primary-bound best_fits.tsv for the wide run's nested rule")
ap.add_argument("--out", required=True)
ap.add_argument("--no_gim", action="store_true")
a = ap.parse_args()
OUT = os.path.join(os.getcwd(), os.environ["FD_OUT"])
spec = importlib.util.spec_from_file_location("fitmod", os.path.join(HERE, "..", "..", "fixed_differences", "07_fit_models.py"))
fitmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fitmod)
ext_models.register(fitmod)
LADDER = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]
MU = (4e-9, 1e-8)
EPS = 1e-4
FIXED = {"pair", "model", "target", "seed", "start", "LL", "theta", "seconds"}

best = {}
for r in csv.DictReader(open(os.path.join(OUT, "fits/best_fits.tsv")), delimiter="\t"):
    best[(r["pair"], r["model"])] = r
ref = {}
if a.refdir:
    for r in csv.DictReader(open(a.refdir), delimiter="\t"):
        ref[(r["pair"], r["model"])] = r


def kv(s):
    return {k: float(x) for k, x in (t.split("=") for t in s.split(";"))} if s else {}


def mutag(mu):
    return f"mu{mu:.0e}".replace("e-0", "e-")


def gim_cov(pair, model, pnames, popt, data, pk):
    """The 08 GIM call (same bootstrap spectra and function), returning cov over (params, theta) and the SEs."""
    func_ex, pts_l, _, _, _, fixed = fitmod.model_setup(model, data.sample_sizes)
    boots = fitmod.bootstrap_spectra(pk, block_mb=1, folded=fitmod.is_folded(model))
    se, gim, hess = dadi.Godambe.GIM_uncert(func_ex, pts_l, boots, popt, data, multinom=True, return_GIM=True)
    return np.linalg.inv(gim), np.asarray(se)


def delta_se(fn, q, cov):
    q = np.asarray(q, float)
    g = np.zeros(len(q))
    for i in range(len(q)):
        h = 1e-6 * abs(q[i]) if q[i] != 0 else 1e-9
        up, dn = q.copy(), q.copy()
        up[i] += h
        dn[i] -= h
        g[i] = (fn(up) - fn(dn)) / (2 * h)
    v = float(g @ cov @ g)
    return float(np.sqrt(v)) if v >= 0 else float("nan")


def one(job):
    pair, model = job
    r = best.get((pair, model))
    if r is None:
        return dict(pair=pair, model=model, note="no fit")
    pk = pickle.load(open(os.path.join(OUT, f"sfs/{pair}.pkl"), "rb"))
    n1, n2 = pk["n"]
    i_idx, j_idx = np.meshgrid(np.arange(n1 + 1), np.arange(n2 + 1), indexing="ij")
    D = i_idx / n1 - j_idx / n2
    # ---- per-model admissibility block, as in the drift test's 09_tail.py ----
    folded = fitmod.is_folded(model)
    v = "_alt" if folded else ""
    obs, L_mask = np.asarray(pk[f"data{v}"]), pk[f"L_mask{v}"]
    data = fitmod.data_of(pk, model)
    func_ex, pts_l, pnames, _, _, _ = fitmod.model_setup(model, data.sample_sizes)
    params = {k: float(x) for k, x in (kv_.split("=") for kv_ in r["params"].split(";"))}
    popt = [params[p] for p in pnames]
    m = func_ex(popt, data.sample_sizes, pts_l)
    theta = dadi.Inference.optimal_sfs_scaling(m.fold() if folded else m, data)
    exp = np.ma.filled(theta * m, 0.0)
    ladder = []
    for thr in LADDER:
        s_ = (np.abs(D) >= thr - EPS) if folded else (D >= thr - EPS)
        ladder.append((thr, float(obs[s_].sum()), float(exp[s_].sum())))
    gate = [(t, O, E) for t, O, E in ladder if E >= 10]
    thr_g, O_g, E_g = gate[-1] if gate else (float("nan"),) * 3
    adequate = bool(gate) and (1 / 3) <= O_g / E_g <= 3
    misid = params.get("p_misid", float("nan"))
    converged = r["converged"] == "True"
    # ---- end of admissibility block ----
    converged_08 = converged
    is_ext = model in ext_models.EXT
    # negative expected cells at the best fit: dadi masks them out of the likelihood (ll_per_bin), so a best fit
    # there could be inflated. Declared rule (before any fit): if such cells carry a data fraction above dadi's own
    # warning cutoff (1e-6), the fit is "numerics_invalid" and not converged.
    mm = np.ma.filled(m, 0.0)
    unmasked = ~np.ma.getmaskarray(data)
    mfold = np.ma.filled(m.fold(), 0.0) if folded else mm
    negc = unmasked & (mfold < 0)
    dat = np.ma.filled(data, 0.0)
    n_neg, d_neg = int(negc.sum()), float(dat[negc].sum())
    numerics_invalid = d_neg / float(dat[unmasked].sum()) > 1e-6
    if numerics_invalid and is_ext:      # baseline rows: reported only, drift-test admissibility unchanged
        converged = False
    row = dict(pair=pair, model=model, family="alternative" if is_ext else "baseline",
               analysis="folded" if folded else "unfolded", n_starts=r["n_starts"], LL=r["LL"],
               LL_gap_top3=r["LL_gap_top3"], c_adj=r["c_adj"], gap_adj=r["gap_adj"],
               converged_strict=r["converged_strict"], converged_08=converged_08, nested_ok_08=r["nested_ok"])
    nested_ok = True
    below = []
    if is_ext:
        f = "f" if folded else ""
        subs = [(s + f, best.get((pair, s + f))) for s in ext_models.SUBMODELS[model.rstrip("f")]]
        if ref:
            subs.append((model + "@primary", ref.get((pair, model))))
        for s, rs in subs:
            if rs is not None and float(r["LL"]) < float(rs["LL"]):
                nested_ok = False
                below.append(f"{s}:{float(rs['LL']) - float(r['LL']):.1f}")
        converged = converged and nested_ok
    adm = converged and adequate and (folded or misid <= 0.05)
    row.update(nested_ok_ext=nested_ok if is_ext else "", below_submodels=";".join(below),
               n_neg_cells=n_neg, data_in_neg_cells=f"{d_neg:.6g}", numerics_invalid=numerics_invalid, converged=converged,
               p_misid=f"{misid:.4f}", gate_threshold=thr_g, gate_O=f"{O_g:.1f}", gate_E=f"{E_g:.1f}",
               gate_O_over_E=f"{O_g / E_g:.3f}" if gate else "nan", adequate=adequate, admissible=adm,
               theta=f"{theta:.1f}", L_mask=f"{L_mask:.0f}", n_BP=n1, n_TP=n2,
               ladder=";".join(f"{t}:{O:.1f}/{E:.1f}" for t, O, E in ladder), params=r["params"], se_08=r.get("se", ""))
    se08 = kv(r.get("se", ""))
    if not is_ext:
        for mu in MU:          # physical conversion as simulation/01_fitted_history.py (baseline rows, context only)
            Na = theta / (4 * mu * L_mask)
            if model.rstrip("f") == "M3":
                nbs, nbe, nts, nte = params["s_BP"] * Na, params["nu_BP"] * Na, (1 - params["s_BP"]) * Na, params["nu_TP"] * Na
            else:
                nbs = nbe = params["nu_BP"] * Na; nts = nte = params["nu_TP"] * Na
            mb = params.get("m_BP_from_TP", params.get("m", 0.0)); mt = params.get("m_TP_from_BP", params.get("m", 0.0))
            tag = mutag(mu)
            row.update({f"N_anc_{tag}": f"{Na:.6g}", f"N_BP_start_{tag}": f"{nbs:.6g}", f"N_BP_end_{tag}": f"{nbe:.6g}",
                        f"N_TP_start_{tag}": f"{nts:.6g}", f"N_TP_end_{tag}": f"{nte:.6g}",
                        f"T_gen_{tag}": f"{params['T'] * 2 * Na:.6g}",
                        f"m_BP_from_TP_pergen_{tag}": f"{mb / (2 * Na):.6g}", f"m_TP_from_BP_pergen_{tag}": f"{mt / (2 * Na):.6g}"})
        return row
    # ---- Stage 2b: best start, its box, bound flags, exact T_c window ----
    seed_rows = []
    for fpath in glob.glob(os.path.join(OUT, f"fits/{pair}.{model}.data.seed*.tsv")):
        for rr in csv.DictReader(open(fpath), delimiter="\t"):
            seed_rows.append(rr)
    bs = max(seed_rows, key=lambda x: float(x["LL"]))
    side = json.load(open(os.path.join(OUT, f"fits/bounds/{pair}.{model}.data.seed{bs['seed']}.json")))
    box = side["box"]
    at = []
    for p in pnames:
        lo, hi = box[p]
        val = params[p]
        if (lo == 0 and val < 1e-3) or (lo > 0 and val <= lo * 1.01):
            at.append(f"{p}=lower")
        if val >= hi / 1.01:
            at.append(f"{p}=upper")
    x = ext_models.tc_x(theta, params["T_c"], L_mask)
    xlo, xhi = ext_models.X_WINDOW[a.bounds]
    if x < xlo / 1.01:
        win = "below_window"
    elif x <= xlo * 1.01:
        win = "at_lower"
    elif x > xhi * 1.01:
        win = "above_window"
    elif x >= xhi / 1.01:
        win = "at_upper"
    else:
        win = "inside"
    if win != "inside" and not any(s.startswith("T_c=") for s in at):
        at.append(f"T_c={win}(exact)")
    row.update(best_seed=bs["seed"], theta_ref=side["theta_ref"], bounds_mode=side["bounds_mode"],
               T_c_box=f"{box['T_c'][0]:.6g}-{box['T_c'][1]:.6g}", T_c_x=f"{x:.6g}", T_c_window=win,
               T_c_gen_mu4e9=f"{x / 4e-9:.6g}", T_c_gen_mu1e8=f"{x / 1e-8:.6g}", at_bound=";".join(at) or "none")
    # ---- SEs (GIM as 08) and derived / physical quantities ----
    cov = None
    # 08 computes the GIM at the best start's full-precision parameters (seed file), not at the 6-significant-digit
    # best_fits values; the same point is used here, so the SEs equal 08's and the delta method is taken there
    popt_full = [float(bs[p]) for p in pnames]
    row["best_start_matches_best_fits"] = all(f"{u:.6g}" == f"{params[p]:.6g}" for u, p in zip(popt_full, pnames))
    m_full = func_ex(popt_full, data.sample_sizes, pts_l)
    theta_full = dadi.Inference.optimal_sfs_scaling(m_full.fold() if folded else m_full, data)
    if not a.no_gim:
        try:
            cov, se = gim_cov(pair, model, pnames, popt_full, data, pk)
            dif = max(abs(se[i] - se08[p]) / max(abs(se08[p]), 1e-300) for i, p in enumerate(pnames) if p in se08) if se08 else float("nan")
            row["se_check_vs_08_maxrel"] = f"{dif:.3g}"
        except Exception as e:
            row["gim_error"] = str(e)
    for i, p in enumerate(pnames):
        row[p] = f"{params[p]:.6g}"
        row[p + "_se"] = f"{se08[p]:.4g}" if p in se08 else ""
    q0 = popt_full + [theta_full]
    ip = {p: i for i, p in enumerate(pnames)}

    def pdict(q):
        return {p: q[i] for p, i in ip.items()}
    dnames = list(ext_models.derived(model, params))
    for dn in dnames:
        fn = (lambda q, dn=dn: ext_models.derived(model, pdict(q))[dn])
        row[dn] = f"{fn(q0):.6g}"
        row[dn + "_se"] = f"{delta_se(fn, q0, cov):.4g}" if cov is not None else ""
    for mu in MU:
        tag = mutag(mu)
        phys0 = ext_models.physical(model, params, theta, L_mask, mu)
        for k in phys0:
            fn = (lambda q, k=k: ext_models.physical(model, pdict(q), q[-1], L_mask, mu)[k])
            row[f"{k}_{tag}"] = f"{phys0[k]:.6g}"
            row[f"{k}_{tag}_se"] = f"{delta_se(fn, q0, cov):.4g}" if cov is not None else ""
    print(pair, model, f"LL {float(r['LL']):.1f} gap_adj {r['gap_adj']} conv {converged} nested_ok {nested_ok} "
          f"misid {misid:.4f} gate {thr_g} O/E {O_g:.1f}/{E_g:.1f} adequate {adequate} admissible {adm} at_bound {row['at_bound']}",
          flush=True)
    return row


jobs = [(p, m) for p in a.pairs.split(",") for m in a.models.split(",")]
n = min(len(jobs), int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
if n > 1:
    with mp.get_context("fork").Pool(n) as pool:
        rows = pool.map(one, jobs, chunksize=1)
else:
    rows = [one(j) for j in jobs]
cols = []
for r in rows:
    for k in r:
        if k not in cols:
            cols.append(k)
with open(a.out + ".tmp", "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in rows:
        o.write("\t".join(str(r.get(k, "")) for k in cols) + "\n")
os.replace(a.out + ".tmp", a.out)
