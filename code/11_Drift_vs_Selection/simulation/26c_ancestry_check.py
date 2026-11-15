#!/usr/bin/env python3
"""26c_ancestry_check.py -- is the 2009->2019 frequency change ancestry rather than drift?
(temporal ancestry check; the unphased projection of the Supplementary Note)

All ten 2019 Apassa fish are wobbles (2 of 20 in 2009), mostly from a site ~0.8 km from the 2009 site.
Reads the genotype tables of 26a_extract_temporal.sh and reuses the loader and site filters of
26b_contemporary_ne.py.

1. Regression: across the sites kept for a population's temporal estimate, the 2009->2019 change in that
   population against the 2009 frequency difference from the other population. If the later sample is
   a mix, (1 - m) own + m other, the slope estimates m; under drift alone the slope is about 0.
2. Per fish: on sites where the 2009 samples of the two populations differ by >= 0.3, a projection score
   q = sum((g/2 - p_APA)(p_BAM - p_APA)) / sum((p_BAM - p_APA)^2), with each fish left out of its own
   group's frequency. About 0 for an Apassa-like genome, about 1 for a Bambomo-like one.

Usage (repo root): 26c_ancestry_check.py <vcf> <26_temporal dir> <out.tsv>
"""
import importlib.util, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("ne26b", os.path.join(HERE, "26b_contemporary_ne.py"))
ne = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ne)


def main():
    vcf, d, out = sys.argv[1:4]
    labels, data = ne.load(vcf, d)
    names = ne.subprocess.run(f"bcftools view -h -S {d}/samples.txt {vcf} 2>/dev/null | grep '^#CHROM'",
                              shell=True, capture_output=True, text=True).stdout.split()[9:]
    G = np.vstack([g for _, g in data.values()])
    idx = {g: np.flatnonzero(labels == g) for g in ("BAM2009", "BAM2019", "APA2009", "APA2019")}
    ok = G >= 0
    a = {g: np.where(ok[:, c], G[:, c], 0).sum(1).astype(float) for g, c in idx.items()}
    n = {g: ok[:, c].sum(1).astype(float) for g, c in idx.items()}
    with np.errstate(invalid="ignore", divide="ignore"):
        p = {g: a[g] / (2 * n[g]) for g in idx}
    lines = []
    for pop, other in (("APA", "BAM"), ("BAM", "APA")):
        c0, c1 = f"{pop}2009", f"{pop}2019"
        with np.errstate(invalid="ignore", divide="ignore"):
            pool = (a[c0] + a[c1]) / (2 * (n[c0] + n[c1]))
        keep = ((n[c0] >= ne.CALL_T * idx[c0].size) & (n[c1] >= ne.CALL_T * idx[c1].size)
                & (n[f"{other}2009"] >= ne.CALL_T * idx[f"{other}2009"].size)
                & (pool >= ne.MAF) & (pool <= 1 - ne.MAF))
        dd = (p[f"{other}2009"] - p[c0])[keep]
        dl = (p[c1] - p[c0])[keep]
        slope = float((dd * dl).sum() / (dd * dd).sum())
        r = float(np.corrcoef(dd, dl)[0, 1])
        lines.append(f"regression\t{pop}\ttoward_{other}\tsites={keep.sum()}\tslope_m={slope:.4f}\tpearson_r={r:.4f}")
        print(lines[-1], flush=True)
    # per-fish projection score on informative sites
    base_ok = ((n["APA2009"] >= 0.8 * idx["APA2009"].size) & (n["BAM2009"] >= 0.8 * idx["BAM2009"].size))
    with np.errstate(invalid="ignore"):
        inf_sites = base_ok & (np.abs(p["BAM2009"] - p["APA2009"]) >= 0.3)
    Gi, oki = G[inf_sites], ok[inf_sites]
    lines.append(f"informative_sites\t{inf_sites.sum()}")
    for g, cols in idx.items():
        for j in cols:
            pa, pb = p["APA2009"][inf_sites].copy(), p["BAM2009"][inf_sites].copy()
            if g in ("APA2009", "BAM2009"):             # leave the fish out of its own group's frequency
                aa, nn = a[g][inf_sites], n[g][inf_sites]
                gj = np.where(oki[:, j], Gi[:, j], 0); cj = oki[:, j].astype(float)
                loo = (aa - gj) / (2 * (nn - cj))
                if g == "APA2009": pa = loo
                else: pb = loo
            m = oki[:, j]
            dif = (pb - pa)[m]
            q = float(((Gi[m, j] / 2 - pa[m]) * dif).sum() / (dif * dif).sum())
            lines.append(f"fish\t{g}\t{names[j]}\t{q:.3f}")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    for g in idx:
        qs = [float(l.split("\t")[3]) for l in lines if l.startswith(f"fish\t{g}\t")]
        print(f"{g}: n={len(qs)} q median {np.median(qs):.3f} range {min(qs):.3f} to {max(qs):.3f}", flush=True)


if __name__ == "__main__":
    main()
