#!/usr/bin/env python3
"""write_admissible.py -- hand-off table for Stage 3: one row per Stage 2b alternative
(CB, CBT, CBTM, CBf, CBTf, CBTMf) on BAM x TP1 from gate_ext.tsv, with the admissibility decision and its parts,
LL, Godambe gap, every parameter in dadi units (fitted and derived, with Godambe SEs), physical units for
mu = 1e-8 and 4e-9 (N_anc and all sizes, times, migration rates per generation, with delta-method SEs), p_misid and the
at-bound flags. Written atomically (.tmp then os.replace); then DONE marker (also .tmp then os.replace).
stdlib only. Usage: python3 write_admissible.py <gate_ext.tsv> <out.tsv> <done_marker>
"""
import csv
import os
import sys

src, out, done = sys.argv[1:4]
rows = [r for r in csv.DictReader(open(src), delimiter="\t") if r.get("family") == "alternative"]
head = ["model", "analysis", "admissible", "converged", "converged_08", "nested_ok_ext", "numerics_invalid", "adequate",
        "gate_threshold", "gate_O", "gate_E", "gate_O_over_E", "LL", "c_adj", "gap_adj", "n_starts", "p_misid", "theta",
        "L_mask", "bounds_mode", "at_bound", "T_c_box", "T_c_x", "T_c_window", "T_c_gen_mu4e9", "T_c_gen_mu1e8", "params"]
dadi_cols = ["nu_BP", "nu_TP", "T_pre", "T_c", "r_BP", "r_TP", "m_BP_from_TP", "m_TP_from_BP", "f_mig", "p_misid",
             "nu_BP2", "nu_TP2", "T_s", "T_m"]
phys = ["N_anc", "N_BP_1", "N_BP_2", "N_TP_1", "N_TP_2", "T_s_gen", "T_c_gen", "T_m_gen", "m_BP_from_TP_pergen",
        "m_TP_from_BP_pergen"]
cols = list(head)
for c in dadi_cols:
    if c != "p_misid":
        cols += [c, c + "_se"]
cols += ["p_misid_se"]
for tag in ("mu1e-8", "mu4e-9"):
    for c in phys:
        cols += [f"{c}_{tag}", f"{c}_{tag}_se"]
with open(out + ".tmp", "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in rows:
        o.write("\t".join(r.get(c, "") for c in cols) + "\n")
os.replace(out + ".tmp", out)
with open(done + ".tmp", "w") as o:
    o.write(f"{len(rows)} rows from {src}\n")
os.replace(done + ".tmp", done)
print(f"wrote {out} ({len(rows)} rows) and {done}")
