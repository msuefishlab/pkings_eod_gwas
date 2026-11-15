#!/usr/bin/env python3
"""alt_expected_sfs.py -- dadi's expected joint spectrum (theta * model, drift-test projection 60 x 60, misid wrapper of
07_fit_models.model_setup for unfolded models, folded for folded models) at the best fit of each listed model, from
the Stage 2b fit root (sfs/A1w.pkl, fits/best_fits.tsv). Run inside dadi_tools.sif.
Usage: FD_OUT=<fit root> alt_expected_sfs.py <stage2b code dir> <models,comma> <out.npz>
"""
import csv, importlib.util, os, pickle, sys
import numpy as np
import dadi

code, models, out = sys.argv[1], sys.argv[2].split(","), sys.argv[3]
sys.path.insert(0, code)
import ext_models  # noqa: E402
spec = importlib.util.spec_from_file_location("fitmod", os.path.join(code, "..", "..", "fixed_differences", "07_fit_models.py"))
fitmod = importlib.util.module_from_spec(spec); spec.loader.exec_module(fitmod)
ext_models.register(fitmod)
OUT = os.environ["FD_OUT"]
best = {(r["pair"], r["model"]): r for r in csv.DictReader(open(os.path.join(OUT, "fits/best_fits.tsv")), delimiter="\t")}
pk = pickle.load(open(os.path.join(OUT, "sfs/A1w.pkl"), "rb"))
res = {}
for model in models:
    r = best[("A1w", model)]
    folded = fitmod.is_folded(model)
    v = "_alt" if folded else ""
    data = fitmod.data_of(pk, model)
    func_ex, pts_l, pnames, _, _, _ = fitmod.model_setup(model, data.sample_sizes)
    params = {k: float(x) for k, x in (kv.split("=") for kv in r["params"].split(";"))}
    m = func_ex([params[p] for p in pnames], data.sample_sizes, pts_l)
    mm = m.fold() if folded else m
    theta = dadi.Inference.optimal_sfs_scaling(mm, data)
    res[f"{model}|exp"] = np.ma.filled(theta * mm, 0.0)
    res[f"{model}|mask"] = np.ma.getmaskarray(data)
    res[f"{model}|L_mask"] = np.array(pk[f"L_mask{v}"])
    res[f"{model}|theta"] = np.array(theta)
    res[f"{model}|obs"] = np.ma.filled(data, 0.0)
    print(model, "theta", theta, "L_mask", pk[f"L_mask{v}"], "folded", folded, flush=True)
np.savez(out, **res)
