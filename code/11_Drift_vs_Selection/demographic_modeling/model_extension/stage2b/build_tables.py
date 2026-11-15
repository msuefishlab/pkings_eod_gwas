#!/usr/bin/env python3
"""build_tables.py -- assemble the Stage 2b deliverable tables from gate_ext.tsv / lrt_ext.tsv (no new statistic).
  stage2b_fits.tsv               BAM x TP1, primary bounds: 6 alternatives + the 8 re-collected drift-test baseline rows (context)
  stage2b_lrt.tsv                Godambe-adjusted LRTs CB vs M0, CBT vs CB, CBTM vs CBT (information only)
  stage2b_bound_sensitivity.tsv  primary vs wide-bound fit of each preferred model with a parameter at a bound
  stage2b_variant_fits.tsv       best alternative(s) on TP1_noAdmix / TP1_no017 / BAM_noKin x TP1 (+ their baseline rows)
Usage: python3 build_tables.py <primary_dir> [<wide_dir>] [<var_dir>] --out <dir>
"""
import argparse
import os

import pandas as pd

SAMPLE = {"A1w": "BAM x TP1", "V2w_noAdm10": "BAM x TP1_noAdmix (37)", "V2w_no017": "BAM x TP1_no017 (40)",
          "V2w_noKin": "BAM_noKin (66) x TP1"}
FRONT = ["sample", "pair", "analysis", "model", "family", "bounds_mode", "n_starts", "LL", "c_adj", "gap_adj",
         "converged_strict", "converged_08", "nested_ok_ext", "below_submodels", "numerics_invalid", "n_neg_cells",
         "converged", "p_misid", "gate_threshold", "gate_O", "gate_E", "gate_O_over_E", "adequate", "admissible",
         "at_bound", "theta", "theta_ref", "L_mask", "T_c_box", "T_c_x", "T_c_window", "T_c_gen_mu4e9", "T_c_gen_mu1e8"]
DROP = {"ladder", "se_08", "best_start_matches_best_fits", "se_check_vs_08_maxrel", "params", "nested_ok_08",
        "n_BP", "n_TP", "data_in_neg_cells", "best_seed"}


def tidy(g, bounds_default="primary"):
    g = g.copy()
    g.insert(0, "sample", g["pair"].map(SAMPLE))
    if "bounds_mode" not in g:
        g["bounds_mode"] = ""
    g.loc[g["family"] == "baseline", "bounds_mode"] = "drift test"
    front = [c for c in FRONT if c in g]
    rest = [c for c in g.columns if c not in front and c not in DROP]
    keep = front + rest + [c for c in ("params", "se_08", "ladder", "best_seed", "se_check_vs_08_maxrel") if c in g]
    return g[keep]


ap = argparse.ArgumentParser()
ap.add_argument("primary")
ap.add_argument("--wide", default=None)
ap.add_argument("--var", default=None)
ap.add_argument("--out", required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
P = pd.read_csv(os.path.join(a.primary, "gate_ext.tsv"), sep="\t")
order = {m: i for i, m in enumerate(["CB", "CBT", "CBTM", "CBf", "CBTf", "CBTMf", "M0", "M1", "M2", "M3",
                                     "M0f", "M1f", "M2f", "M3f"])}
P = P.sort_values("model", key=lambda s: s.map(order))
tidy(P).to_csv(os.path.join(a.out, "stage2b_fits.tsv"), sep="\t", index=False)
L = pd.read_csv(os.path.join(a.primary, "lrt_ext.tsv"), sep="\t")
L.insert(0, "sample", L["pair"].map(SAMPLE))
L["reading"] = "information only; selects nothing; nominal p (nuisance parameter unidentified under H0)"
L.to_csv(os.path.join(a.out, "stage2b_lrt.tsv"), sep="\t", index=False)
if a.wide:
    Wd = pd.read_csv(os.path.join(a.wide, "gate_ext.tsv"), sep="\t")
    pref = P[P["model"].isin(Wd["model"])]
    B = pd.concat([tidy(pref), tidy(Wd)], ignore_index=True)
    B.insert(1, "fit", ["primary bounds"] * len(pref) + ["wide bounds (sensitivity)"] * len(Wd))
    B.to_csv(os.path.join(a.out, "stage2b_bound_sensitivity.tsv"), sep="\t", index=False)
if a.var:
    V = pd.read_csv(os.path.join(a.var, "gate_ext.tsv"), sep="\t")
    V = V.sort_values(["pair", "model"], key=lambda s: s.map(order) if s.name == "model" else s)
    ref = P[P["model"].isin(V.loc[V["family"] == "alternative", "model"])]
    tidy(pd.concat([ref, V], ignore_index=True)).to_csv(os.path.join(a.out, "stage2b_variant_fits.tsv"), sep="\t", index=False)
print("written to", a.out)
