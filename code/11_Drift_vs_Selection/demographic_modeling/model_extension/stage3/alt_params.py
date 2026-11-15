#!/usr/bin/env python3
"""alt_params.py -- convert the Stage 2b admissible alternatives to simulation parameters
(simulation/01_fitted_history.py conventions, extended to the new parameters).
  N_anc = theta / (4 mu L_mask)   (L_mask of the fitted spectrum; folded models: the folded-analysis L_mask)
  sizes = nu * N_anc: BP nu_BP (split -> T_c), r_BP * nu_BP (T_c -> present); TP nu_TP, r_TP * nu_TP (CB: r_TP = 1)
  times: T_s,gen = 2 N_anc (T_pre + T_c); T_c,gen = 2 N_anc T_c; T_m,gen = 2 N_anc f_mig T_s (CBTM)
  m per generation = m / (2 N_anc); dadi m12 = into BP from TP (msprime source BP -> dest TP, as 23b)
The conversion is recomputed here from the dadi-unit parameters and theta and compared with the physical columns of
admissible_alternatives.tsv (relative difference must be <= 1e-4), then written as
  alt_params.tsv  (model, analysis, p_misid, phys_json [mu = 1e-8 arm, used by the simulation])
  alt_params.json ({"<model>|1e-08": {...}, "<model>|4e-09": {...}} for the recent-Ne diagnostic)
Usage: alt_params.py <admissible_alternatives.tsv> <out_prefix>
"""
import csv, json, sys

src, outp = sys.argv[1:3]
PH = ["N_anc", "N_BP_1", "N_BP_2", "N_TP_1", "N_TP_2", "T_s_gen", "T_c_gen", "T_m_gen", "m_BP_from_TP_pergen",
      "m_TP_from_BP_pergen"]
rows, js, report = [], {}, []
for r in csv.DictReader(open(src), delimiter="\t"):
    if r["admissible"] != "True":
        report.append(f"{r['model']}: not admissible (skipped)")
        continue
    p = {k: float(v) for k, v in (kv.split("=") for kv in r["params"].split(";") if kv)}
    theta, L = float(r["theta"]), float(r["L_mask"])
    base = r["model"].rstrip("f")
    for mu, tag in ((1e-8, "mu1e-8"), (4e-9, "mu4e-9")):
        Na = theta / (4 * mu * L)
        Ts = p["T_pre"] + p["T_c"]
        q = dict(N_anc=Na, N_BP_1=p["nu_BP"] * Na, N_BP_2=p["r_BP"] * p["nu_BP"] * Na, N_TP_1=p["nu_TP"] * Na,
                 N_TP_2=p.get("r_TP", 1.0) * p["nu_TP"] * Na, T_s_gen=2 * Na * Ts, T_c_gen=2 * Na * p["T_c"])
        if base == "CBTM":
            q.update(T_m_gen=2 * Na * p["f_mig"] * Ts, m_BP_from_TP_pergen=p["m_BP_from_TP"] / (2 * Na),
                     m_TP_from_BP_pergen=p["m_TP_from_BP"] / (2 * Na))
        else:
            q.update(m_BP_from_TP_pergen=0.0, m_TP_from_BP_pergen=0.0)
        worst = 0.0
        for k, v in q.items():
            t = r.get(f"{k}_{tag}", "")
            if t not in ("", "nan"):
                worst = max(worst, abs(float(t) - v) / max(abs(v), 1e-300))
        report.append(f"{r['model']} {tag}: max rel diff vs table {worst:.2e}")
        assert worst <= 1e-4, f"{r['model']} {tag}: conversion disagrees with admissible_alternatives.tsv ({worst})"
        js[f"{r['model']}|{mu:g}"] = q
    pm = r["p_misid"] if r["analysis"] == "unfolded" else "0"
    rows.append(dict(model=r["model"], analysis=r["analysis"], p_misid=pm, phys_json=json.dumps(js[f"{r['model']}|1e-08"]),
                     at_bound=r.get("at_bound", ""), T_c_window=r.get("T_c_window", "")))
with open(outp + ".tsv", "w") as fo:
    w = csv.DictWriter(fo, fieldnames=["model", "analysis", "p_misid", "phys_json", "at_bound", "T_c_window"],
                       delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
json.dump(js, open(outp + ".json", "w"), indent=1)
print("\n".join(report))
print("admissible alternatives:", ",".join(r["model"] for r in rows) or "none")
