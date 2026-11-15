#!/usr/bin/env python3
"""07_match.py -- simulation-calibrated reading of the real fish.

Matches every real fish to the simulated cells (05_simulate.py, sim/fish.tsv) with la_lib.match:
  phased   (q_m, m_hat, log10 T_hat, F_long) from hap_primary/fits.tsv
  diploid  (q_m, m_hat_dip, log10 T_hat_dip) from dip_phasedgeno/ and dip_unphased/fits.tsv
Matched T = best cell's T; interval = range of T over compatible cells (|z| <= 3 on every statistic).
Coarse class: recent (T <= 50), older (100-3000), no resolvable tract structure (inf or no
compatible cell), no admixture (m = 0 cell); "between" when the compatible cells span classes. Final class:
the phased match, unless the diploid (phased-genotype) matched T is > 2-fold below the phased one for the
median 2019 Apassa fish (switch-error rule), then the diploid match.

Also: Pearson r of q_m (and m_hat) with the 26c score; the site check (2019-017 vs 2019-018/021; exact
two-sided Mann-Whitney on q_m, matched m, matched T); panel sensitivity (does every 2019 Apassa fish exceed every
2009 Apassa TP fish in q_m, m_hat and F_long under the primary, full and BEN panels); tract sharing among the
2019 Apassa fish and, for reference, among the 2009 Apassa fish outside the panels.

Writes <out_dir>/results/: per_fish.tsv, site_check.tsv, panel_sensitivity.tsv, tract_sharing.tsv, summary.txt.

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/07_match.py <out_dir>
"""
import csv, os, sys
import numpy as np
from scipy.stats import mannwhitneyu, pearsonr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import la_lib as L


def read(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def tnum(x):
    return float(x)


def cclass(T, m):
    T, m = float(T), float(m)
    if m == 0:
        return "no_admixture"
    if np.isinf(T):
        return "no_tract_structure"
    return "recent" if T <= 50 else "older"


def classify(best, compat):
    if not compat:
        return "no_tract_structure(no_cell_fits)"
    cl = sorted({cclass(*c) for c in compat})
    return cl[0] if len(cl) == 1 else "between:" + "/".join(cl)


def sim_rows(out):
    rows = read(f"{out}/sim/fish.tsv")
    for r in rows:
        r["T"] = tnum(r["T"]); r["m"] = float(r["m"])
    return rows


def main():
    out = sys.argv[1]
    od = f"{out}/results"; os.makedirs(od, exist_ok=True)
    srows = sim_rows(out)
    cells = {k: L.cell_stats(srows, keys) for k, keys in L.MATCH_KEYS.items()}
    hap = {r["sample"]: r for r in read(f"{out}/hap_primary/fits.tsv")}
    dips = {src: {r["sample"]: r for r in read(f"{out}/dip_{src}/fits.tsv")} for src in ("phasedgeno", "unphased")}
    res = []
    for s, h in hap.items():
        x = dict(q_m=float(h["q_m"]), m_hat=float(h["m"]), logT=float(h["logT"]), F_long=float(h["F_long"]))
        best, compat, _ = L.match(x, cells["phased"], L.MATCH_KEYS["phased"])
        r = dict(sample=s, group=h["group"], site=h["site"], q26c=h["q26c"], q_m=x["q_m"], m_hat=x["m_hat"],
                 T_hat=10 ** x["logT"], F_long=x["F_long"], n_tracts=int(h["n_tracts"]),
                 median_tract=float(h["median_tract_bp"]), longest=float(h["longest_tract_bp"]),
                 frac_tract=float(h["frac_tract"]),
                 match_T=best[0], match_m=best[1],
                 match_T_lo=min((c[0] for c in compat), default=np.nan),
                 match_T_hi=max((c[0] for c in compat), default=np.nan),
                 match_m_lo=min((c[1] for c in compat), default=np.nan),
                 match_m_hi=max((c[1] for c in compat), default=np.nan),
                 n_compat=len(compat), class_phased=classify(best, compat))
        for src, D in dips.items():
            d = D[s]
            y = dict(q_m=float(d["q_m"]), m_hat_dip=float(d["m"]), logT_dip=float(d["logT"]))
            b, cp, _ = L.match(y, cells["diploid"], L.MATCH_KEYS["diploid"])
            r[f"q_m_{src}"] = y["q_m"]; r[f"m_hat_{src}"] = y["m_hat_dip"]; r[f"T_hat_{src}"] = 10 ** y["logT_dip"]
            r[f"match_T_{src}"] = b[0]; r[f"match_m_{src}"] = b[1]
            r[f"match_T_lo_{src}"] = min((c[0] for c in cp), default=np.nan)
            r[f"match_T_hi_{src}"] = max((c[0] for c in cp), default=np.nan)
            r[f"class_{src}"] = classify(b, cp)
            if src == "unphased":
                r["call_rate_unphased"] = float(d["call_rate"]); r["concordance"] = float(d["concordance_with_phased"])
        res.append(r)
    a19 = [r for r in res if r["group"] == "APA2019_Wobble"]
    ph = np.median([r["match_T"] for r in a19]); dp = np.median([r["match_T_phasedgeno"] for r in a19])
    use_dip = bool(dp * 2 < ph)
    for r in res:
        r["class_final"] = r["class_phasedgeno"] if use_dip else r["class_phased"]
        r["T_r1.0"] = r["match_T"] * 1.4 / 1.0; r["T_r3.0"] = r["match_T"] * 1.4 / 3.0
    order = ["APA2019_Wobble", "APA2009_Wobble", "APA2009_TP", "BEN2023_TP", "BAM2009_BP", "BAM2009_Wobble",
             "BAM2009_TP", "BAM2019_Wobble"]
    res.sort(key=lambda r: (order.index(r["group"]), r["site"], r["sample"]))
    keys = list(res[0])
    with open(f"{od}/per_fish.tsv", "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in res:
            fh.write("\t".join(f"{r[k]:.4g}" if isinstance(r[k], (float, np.floating)) else str(r[k]) for k in keys) + "\n")

    lines = [f"switch-error rule: median matched T phased {ph:g}, diploid(phased genotypes) {dp:g}; "
             f"use diploid = {use_dip}"]
    q26 = np.array([float(r["q26c"]) for r in res]); qm = np.array([r["q_m"] for r in res])
    ok = np.isfinite(q26)
    lines.append(f"Pearson q_m vs q26c, {ok.sum()} fish with a 26c score: r = {pearsonr(qm[ok], q26[ok])[0]:.4f}")
    apa = [r for r in res if r["group"].startswith("APA")]
    lines.append(f"Pearson q_m vs q26c, {len(apa)} Apassa fish: r = "
                 f"{pearsonr([r['q_m'] for r in apa], [float(r['q26c']) for r in apa])[0]:.4f}")

    # site check
    g7 = [r for r in a19 if r["site"] == "2019-017"]; g3 = [r for r in a19 if r["site"] != "2019-017"]
    sc = ["statistic\tn_017\tn_018_021\tvalues_017\tvalues_018_021\tp_exact_two_sided\toverlap"]
    for k in ("q_m", "match_m", "match_T", "m_hat", "F_long"):
        u, v = [r[k] for r in g7], [r[k] for r in g3]
        p = mannwhitneyu(u, v, alternative="two-sided", method="exact").pvalue
        ov = not (min(u) > max(v) or max(u) < min(v))
        sc.append(f"{k}\t{len(u)}\t{len(v)}\t{','.join(f'{x:.4g}' for x in u)}\t{','.join(f'{x:.4g}' for x in v)}\t"
                  f"{p:.4f}\t{int(ov)}")
    open(f"{od}/site_check.tsv", "w").write("\n".join(sc) + "\n")

    # panel sensitivity
    ps = ["panel\tstatistic\tAPA2019_min\tAPA2019_median\tAPA2009_TP_max\tAPA2009_TP_median\tn_2009_TP\t"
          "all_2019_above_all_2009TP"]
    for panel in ("primary", "full", "ben"):
        F = read(f"{out}/hap_{panel}/fits.tsv")
        x19 = [r for r in F if r["group"] == "APA2019_Wobble"]; x09 = [r for r in F if r["group"] == "APA2009_TP"]
        for k, col in (("q_m", "q_m"), ("m_hat", "m"), ("F_long", "F_long")):
            u = [float(r[col]) for r in x19]; v = [float(r[col]) for r in x09]
            ps.append(f"{panel}\t{k}\t{min(u):.4f}\t{np.median(u):.4f}\t{max(v):.4f}\t{np.median(v):.4f}\t{len(v)}\t"
                      f"{int(min(u) > max(v))}")
    open(f"{od}/panel_sensitivity.tsv", "w").write("\n".join(ps) + "\n")

    # tract sharing (fish-level coverage on 10-kb bins)
    tr = read(f"{out}/hap_primary/tracts.tsv")
    M = L.load_markers(out, "primary"); gw = M["gw"].astype(bool)
    seg = L.segments(M["chrom"][gw], M["pos"][gw], L.excluded_regions())
    spans = [(int(M["chrom"][gw][seg == s][0]), M["pos"][gw][seg == s].min(), M["pos"][gw][seg == s].max())
             for s in np.unique(seg)]
    BIN = 10_000
    binsets = {}
    for c, a, b in spans:
        for x in range(int(a) // BIN, int(b) // BIN + 1):
            binsets[(c, x)] = len(binsets)
    nb = len(binsets)
    sh = ["set\tn_fish\tobserved_shared_fraction\texpected_if_independent\tmean_fish_coverage"]
    sets = {"APA2019_Wobble": [r["sample"] for r in a19],
            "APA2009_outside_panels": [r["sample"] for r in res if r["group"] in ("APA2009_TP", "APA2009_Wobble")]}
    for name, fs in sets.items():
        cov = np.zeros((len(fs), nb), bool)
        for i, s in enumerate(fs):
            for t in tr:
                if t["sample"] != s:
                    continue
                c = int(t["chrom"][3:])
                for x in range(int(float(t["start"])) // BIN, int(float(t["end"])) // BIN + 1):
                    j = binsets.get((c, x))
                    if j is not None:
                        cov[i, j] = True
        f = cov.mean(1)
        n_other = cov.sum(0)[None, :] - cov
        obs = (cov & (n_other >= 1)).sum() / cov.sum()
        exp = np.sum([cov[i].sum() * (1 - np.prod(np.delete(1 - f, i))) for i in range(len(fs))]) / cov.sum()
        sh.append(f"{name}\t{len(fs)}\t{obs:.4f}\t{exp:.4f}\t{f.mean():.4f}")
    open(f"{od}/tract_sharing.tsv", "w").write("\n".join(sh) + "\n")
    open(f"{od}/summary.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
