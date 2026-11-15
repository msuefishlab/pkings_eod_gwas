#!/usr/bin/env python3
"""write_sim_index.py -- model and ratio indices, written before any Stage 3 simulation.
Seed rule: genome / chunk replicate g of model x at ratio index k: seed = 20261109 + 1,000,000 i_model + 10,000 i_ratio + g.
Usage: write_sim_index.py <params.as_run.agree.json> <agree best_fits.tsv> <out.tsv>
"""
import csv, json, sys

pj, bf, out = sys.argv[1:4]
P = json.load(open(pj))
pm = {}
for r in csv.DictReader(open(bf), delimiter="\t"):
    if r["pair"] == "A1w":
        kv = dict(x.split("=") for x in r["params"].split(";") if x)
        pm[r["model"]] = kv.get("p_misid", "")
rows = []
base = [("M0", "unfolded"), ("M2", "unfolded"), ("M3", "unfolded"), ("M0f", "folded"), ("M1f", "folded"),
        ("M2f", "folded"), ("M3f", "folded")]
for i, (m, an) in enumerate(base, 1):
    assert f"{m}|1e-08" in P
    p = pm[m] if an == "unfolded" else "0"
    rows.append(dict(kind="model", model=m, i_model=i, analysis=an, source="baseline", p_misid=p,
                     phys_json=json.dumps(P[f"{m}|1e-08"]), status="baseline (drift-test agreement run, A1w, mu = 1e-8)",
                     value="", use=""))
for i, (m, an) in enumerate([("CB", "unfolded"), ("CBT", "unfolded"), ("CBTM", "unfolded"), ("CBf", "folded"),
                             ("CBTf", "folded"), ("CBTMf", "folded")], 8):
    rows.append(dict(kind="model", model=m, i_model=i, analysis=an, source="alternative", p_misid="", phys_json="",
                     status="reserved; used only if admissible in stage2b/admissible_alternatives.tsv (parameters then "
                            "in stage3/alt_params.tsv)", value="", use=""))
GRID = [0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 7.5, 10]
for k, v in enumerate(GRID, 1):
    rows.append(dict(kind="ratio", model="", i_model="", analysis="", source="", p_misid="", phys_json="", status="",
                     value=v, use=f"i_ratio {k}: simulation/23e grid value (physical r/mu): reference families (g 1-100 "
                                  f"calibration chunks, g 101-200 all-M-sum chunks) and held-out scoring family "
                                  f"(g 1001-1005 held-out genome sets)"))
rows.append(dict(kind="ratio", model="", i_model="", analysis="", source="", p_misid="", phys_json="", status="",
                 value="fitted", use="i_ratio 20: the model's BAM calibration-fitted r/mu (primary row), fresh "
                                     "whole M-sum replicate genomes g 1-100 (P1 held-out, P3, supporting diagnostics)"))
rows.append(dict(kind="ratio", model="", i_model="", analysis="", source="", p_misid="", phys_json="", status="",
                 value="truth", use="i_ratio 30: false-rejection-rate calibration truth r/mu (calibration-fitted, or 1.4 at a grid edge), "
                                    "pseudo-observed whole M-sum genomes g 1-100"))
rows.append(dict(kind="ratio", model="", i_model="", analysis="", source="", p_misid="", phys_json="", status="",
                 value="1.4", use="i_ratio 99: step-1 pipeline test only (one M2 genome, g 1; timing and format check); "
                                  "enters no statistic"))
rows.append(dict(kind="note", model="", i_model="", analysis="", source="", p_misid="", phys_json="", status="",
                 value="", use="per block j (0-based row of mask_blocks.tsv): SeedSequence([seed, j]) -> ancestry, "
                               "mutation seeds; default_rng([seed, j, k]) k = 1 misid, 2 gqdrop, 3 eps0.005, 4 eps0.01. "
                               "Pilot reproduction uses the pilot's own seeds."))
cols = ["kind", "model", "i_model", "analysis", "source", "p_misid", "value", "use", "status", "phys_json"]
with open(out + ".tmp", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader(); w.writerows(rows)
import os  # noqa: E402
os.replace(out + ".tmp", out)
print("wrote", out, len(rows), "rows; p_misid", {m: pm.get(m) for m, a in base if a == "unfolded"})
