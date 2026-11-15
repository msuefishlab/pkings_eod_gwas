#!/usr/bin/env python3
"""23d_score_clustering.py -- score a clustering arm of the LD-decay test the way table 15
(15_clustering_matched_window_agree.tsv, the matched-window clustering null) was scored.

For each history directory of an arm (04c_worker_msprime_agree.py outputs: chunks/, pos/, params.json): cut every
2-Mb chunk into non-overlapping windows of 353,902 bp (Peaks 1 and 4; sites at delta >= delta_peak4) or 287,163 bp
(Peak 5; sites at delta >= delta_peak5), floor(2 Mb / window) windows per chunk, the remainder unused. A window
counts if it holds at least the observed number of sites (table 15: Peak 1 16, Peak 4 11, Peak 5 9).
E = P(window >= observed) x L_scan / window size.

With --compare, each row is also tested against table 15 (15_clustering_matched_window_agree.tsv) at a given mu
arm: two-sided Fisher exact test of the two window counts.

Usage: 23d_score_clustering.py <arm_dir> <out.tsv> [--compare <table15.tsv> <mu>]
"""
import csv, glob, json, math, os, sys
from collections import defaultdict
from scipy.stats import fisher_exact

PEAKS = [(1, 16, 353_902, "ge4"), (4, 11, 353_902, "ge4"), (5, 9, 287_163, "ge5")]
MU = "1e-08"


def score(d):
    par = json.load(open(os.path.join(d, "params.json")))
    model = os.path.basename(d.rstrip("/"))
    chunk_len = int(par["chunk_len"])
    L_scan = par[f"{model}|{MU}"]["L_scan"]
    chunks = set()
    for f in glob.glob(os.path.join(d, "chunks", "chunks_*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            chunks.add((r["task"], int(r["chunk"])))
    counts = {w: defaultdict(lambda: [0, 0]) for w in (353_902, 287_163)}   # (task, chunk, win) -> [ge4, ge5]
    for f in glob.glob(os.path.join(d, "pos", "pos_*.tsv")):
        for r in csv.DictReader(open(f), delimiter="\t"):
            key = (r["task"], int(r["chunk"]))
            if key not in chunks:
                continue
            p = float(r["pos"])
            for w in counts:
                k = int(p // w)
                if k < chunk_len // w:
                    c = counts[w][key + (k,)]
                    c[1] += 1                       # every pos row has delta >= delta_peak5
                    c[0] += int(r["ge_peak4"])
    rows = []
    for peak, obs, w, col in PEAKS:
        n_win = len(chunks) * (chunk_len // w)
        j = 0 if col == "ge4" else 1
        n_ge = sum(1 for c in counts[w].values() if c[j] >= obs)
        mx = max((c[j] for c in counts[w].values()), default=0)
        p = n_ge / n_win if n_win else float("nan")
        rows.append(dict(model=model, r_rate=par["r_rate"], r_over_mu=par["r_rate"] / float(MU), peak=peak,
                         obs_n=obs, window_bp=w, n_chunks=len(chunks), n_sim_windows=n_win,
                         n_sim_windows_ge_obs=n_ge, p_window_ge_obs=p, E_windows_per_genome=p * L_scan / w,
                         sim_window_max=mx))
    return rows


def main():
    arm, out = sys.argv[1], sys.argv[2]
    ref = None
    if "--compare" in sys.argv:
        i = sys.argv.index("--compare")
        tab, mu = sys.argv[i + 1], float(sys.argv[i + 2])
        ref = {}
        for r in csv.DictReader(open(tab), delimiter="\t"):
            if math.isclose(float(r["mu"]), mu):
                ref[(r["model"], int(r["peak"]))] = (int(r["n_sim_windows_ge_obs"]), int(r["n_sim_windows"]),
                                                     float(r["E_windows_per_genome"]))
    rows = []
    for d in sorted(glob.glob(os.path.join(arm, "M*"))):
        rows += score(d)
    for r in rows:
        if ref is not None and (r["model"], r["peak"]) in ref:
            a, n, e = ref[(r["model"], r["peak"])]
            _, pf = fisher_exact([[r["n_sim_windows_ge_obs"], r["n_sim_windows"] - r["n_sim_windows_ge_obs"]],
                                  [a, n - a]])
            r.update(table15_ge=a, table15_windows=n, table15_E=e, fisher_p=pf, control_pass=pf >= 0.001)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        extra = f"  table15 E {r['table15_E']:.3g}  Fisher p {r['fisher_p']:.3g}" if "fisher_p" in r else ""
        print(f"{r['model']:4s} P{r['peak']}  {r['n_sim_windows_ge_obs']:6d}/{r['n_sim_windows']:7d}  "
              f"E {r['E_windows_per_genome']:.3g}{extra}")


if __name__ == "__main__":
    main()
