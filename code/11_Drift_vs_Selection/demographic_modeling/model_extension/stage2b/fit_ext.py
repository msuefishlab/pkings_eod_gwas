#!/usr/bin/env python3
"""fit_ext.py -- run the drift test's 07_fit_models.py (../../fixed_differences/, unchanged) with the Stage 2b
models registered (ext_models.register). Takes every 07 argument plus:
  --bounds primary|wide   primary bounds or the widened-bound sensitivity
  --theta_ref THETA       reference theta that maps the T_c generation window onto a box in dadi units
                          (round 1: the pair's best M0 / M0f theta; refine rounds: theta of the start being refined;
                          nested passes: theta of the embedded sub-model fit)
Writes a sidecar fits/bounds/<pair><suffix>.<model>.data.seed<seed>.json with the box used.
Usage (inside dadi_tools.sif, cwd = run root, FD_OUT = fits root):
  python3 fit_ext.py --pair A1w --model CB --starts 5 --seed 1 --maxtime 2400 --bounds primary --theta_ref 1016700
"""
import argparse
import importlib.util
import json
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ext_models  # noqa: E402

spec = importlib.util.spec_from_file_location("fitmod", os.path.join(HERE, "..", "..", "fixed_differences", "07_fit_models.py"))
fitmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fitmod)

pre = argparse.ArgumentParser(add_help=False)
pre.add_argument("--bounds", default="primary", choices=sorted(ext_models.X_WINDOW))
pre.add_argument("--theta_ref", type=float, required=True)
pre.add_argument("--pair", required=True)
pre.add_argument("--model", required=True)
pre.add_argument("--seed", type=int, default=1)
pre.add_argument("--suffix", default="")
ext, _ = pre.parse_known_args()
if ext.model not in ext_models.EXT:
    sys.exit(f"fit_ext.py fits only the Stage 2b models {sorted(ext_models.EXT)}; use 07_fit_models.py for {ext.model}")
with open(os.path.join(fitmod.OUT, f"sfs/{ext.pair}{ext.suffix}.pkl"), "rb") as fh:
    pk = pickle.load(fh)
L = float(pk["L_mask_alt"] if fitmod.is_folded(ext.model) else pk["L_mask"])
lo, hi = ext_models.register(fitmod, ext.bounds, ext.theta_ref, L)
pnames = ext_models.PNAMES[ext.model.rstrip("f")] + ([] if fitmod.is_folded(ext.model) else ["p_misid"])
os.makedirs(os.path.join(fitmod.OUT, "fits", "bounds"), exist_ok=True)
side = os.path.join(fitmod.OUT, "fits", "bounds", f"{ext.pair}{ext.suffix}.{ext.model}.data.seed{ext.seed}.json")
with open(side + ".tmp", "w") as o:
    json.dump({"pair": ext.pair, "model": ext.model, "seed": ext.seed, "bounds_mode": ext.bounds,
               "theta_ref": ext.theta_ref, "L": L, "x_window": ext_models.X_WINDOW[ext.bounds],
               "box": {p: list(fitmod.BOUNDS[p]) for p in pnames},
               "start": {p: fitmod.START[p] for p in pnames}}, o, indent=1)
os.replace(side + ".tmp", side)
print(f"fit_ext: {ext.model} bounds={ext.bounds} theta_ref={ext.theta_ref:.6g} L={L:.0f} T_c box [{lo:.6g}, {hi:.6g}]",
      flush=True)
# hand every argument except the two Stage 2b ones to the unchanged 07 main()
argv, skip = [sys.argv[0]], False
for i, a in enumerate(sys.argv[1:]):
    if skip:
        skip = False
        continue
    if a in ("--bounds", "--theta_ref"):
        skip = True
        continue
    if a.startswith("--bounds=") or a.startswith("--theta_ref="):
        continue
    argv.append(a)
sys.argv = argv
fitmod.main()
