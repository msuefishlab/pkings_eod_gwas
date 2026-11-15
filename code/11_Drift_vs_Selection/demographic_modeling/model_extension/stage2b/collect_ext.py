#!/usr/bin/env python3
"""collect_ext.py -- run the drift test's 08_collect_fits.py (../../fixed_differences/, which loads 07_fit_models.py
beside it) unchanged with the Stage 2b models registered. Two run-time additions, no edits to the file:
  1. ext_models.register(mod08.fitmod): model_setup then knows CB, CBT, CBTM (+ f); bounds are not used in collection.
  2. 08 globs fits/*<suffix>.M*.data.seed*.tsv (every drift-test model name starts with M); the glob seen by 08 is widened
     to also return fits/*<suffix>.C*.data.seed*.tsv, appended after 08's own (native-order) result. Nothing else in 08 changes: its Godambe-adjusted convergence,
     Godambe SEs (100 x 1-Mb block bootstrap), its nested_ok rule (M1, M2 only) and its LRTs (M0 vs M1, M2).
     The Stage 2b nested rule and LRTs are applied afterwards by gate_ext.py and lrt_ext.py.
Usage: FD_OUT=<fits root> python3 collect_ext.py --godambe
"""
import glob as _glob
import importlib.util
import os
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ext_models  # noqa: E402

spec = importlib.util.spec_from_file_location("collect08", os.path.join(HERE, "..", "..", "fixed_differences", "08_collect_fits.py"))
mod08 = importlib.util.module_from_spec(spec)
sys.modules["collect08"] = mod08          # so multiprocessing can pickle collect_one by reference
spec.loader.exec_module(mod08)
ext_models.register(mod08.fitmod)
_g = _glob.glob


def _glob_ext(pattern):
    # 08's own glob result first, in its native order (08 keeps the first of LL-tied starts, so the order decides
    # which tied start is reported; the drift-test rows must come out exactly as in its best_fits.tsv), then the Stage 2b files
    found = list(_g(pattern))
    if ".M*." in pattern:
        seen = set(found)
        found += [p for p in _g(pattern.replace(".M*.", ".C*.")) if p not in seen]
    return found


mod08.glob = types.SimpleNamespace(glob=_glob_ext)
mod08.main()
