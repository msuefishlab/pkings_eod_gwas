#!/usr/bin/env python3
"""23f_dense_windows.py -- what are the genome windows as dense as the Bambomo peaks? (LD-decay test, step 2)

Windows: the non-overlapping 354-kb windows with >= 11 sites at delta >= delta_peak4 and the 287-kb windows with
>= 9 sites at delta >= delta_peak5 in 17_observed_window_counts_agree.tsv (the thresholds of table 15), plus the
extended intervals of Peaks 1, 4 and 5 (input_data/10_Introgression/peaks.bed) for reference. For each window, from the unphased VCF and the main
GWAS (PKINGS_ALL_WOB_EXCLUDED, LRT):
  - n_diff: SNPs where the Bambomo (67 A1w fish) and APA + BEN (47) allele frequencies differ by >= 0.78;
  - how many of those were in the GWAS (input filter MAF > 0.05) and passed peak calling (af > 0.15,
    gwas_report.Rmd), and the best -log10 p among them;
  - the Biroundou frequency of the Bambomo-high allele at those SNPs (mean; share >= 0.25);
  - local LD in Bambomo: mean r^2 over 10-50 kb (layout units; ld23_lib) and its percentile among the 200 random
    blocks of 23a_observed_ld_BAM.tsv (high percentile = slow decay = a low-recombination candidate);
  - callable fraction of the block holding the window centre (sfs/A1w.block_scan.tsv).

Usage (repo root): python3 code/11_Drift_vs_Selection/simulation/23f_dense_windows.py <vcf> <out.tsv>
"""
import csv, importlib.util, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ld23_lib as L
spec = importlib.util.spec_from_file_location("obs23a", os.path.join(HERE, "23a_observed_ld_decay.py"))
obs23a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(obs23a)

SIM = "output_data/11_Drift_vs_Selection/simulation"
BLOCKS = "output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree/sfs/A1w.block_scan.tsv"
GWAS = "output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED.lmm2.assoc.txt"
PEAKS = "input_data/10_Introgression/peaks.bed"
POPS = "input_data/05_PopGen/popsfile.txt"
D_MIN = 0.78
LD_BINS = slice(5, 25)          # 10-50 kb
WIN = {"w354k": (353_902, "n_peak4", 11), "w287k": (287_163, "n_peak5", 9)}


def windows():
    out = []
    for r in csv.DictReader(open(os.path.join(SIM, "17_observed_window_counts_agree.tsv")), delimiter="\t"):
        w, col, thr = WIN[r["window_set"]]
        if int(r[col]) >= thr:
            s = int(r["window_index"]) * w
            out.append(dict(window=f"{r['window_set']}:{r['chrom']}:{s}-{s + w}", chrom=r["chrom"], start=s,
                            end=s + w, set=r["window_set"], count=int(r[col]), n_sites=int(r["n_sites"])))
    for line in open(PEAKS):
        p = line.split()
        if p[3] in ("1", "4", "5"):
            out.append(dict(window=f"peak{p[3]}:chr{p[0]}:{p[1]}-{p[2]}", chrom=f"chr{p[0]}", start=int(p[1]),
                            end=int(p[2]), set="peak", count="", n_sites=""))
    return out


def main():
    vcf, out = sys.argv[1], sys.argv[2]
    tmp = os.path.dirname(os.path.abspath(out))
    groups = {}
    for line in open(POPS):
        s, pop = line.split()[:2]
        g = "TP1" if pop in ("APA", "BEN") else pop
        groups.setdefault(g, []).append(s)
    files = {}
    for g in ("BAM", "TP1", "BIR"):
        files[g] = os.path.join(tmp, f"23f_samples_{g}.txt")
        open(files[g], "w").write("\n".join(groups[g]) + "\n")
    frac = {}
    with open(BLOCKS) as fh:
        next(fh)
        for line in fh:
            ch, b, ls = line.split("\t")[:3]
            frac[(ch, int(b))] = int(ls) / 1e6
    # genome distribution of the 10-50 kb mean r^2 over the 200 random blocks
    seg = {}
    for r in csv.DictReader(open(os.path.join(SIM, "23a_observed_ld_BAM.tsv")), delimiter="\t"):
        k = (r["chrom"], r["block"])
        seg.setdefault(k, [np.zeros(L.NBIN), np.zeros(L.NBIN)])
        seg[k][0][int(r["bin"])] += int(r["n_pairs"])
        seg[k][1][int(r["bin"])] += float(r["sum_r2"])
    genome = np.array([t[LD_BINS].sum() / n[LD_BINS].sum() for n, t in seg.values() if n[LD_BINS].sum() > 0])
    wins = windows()
    want = {}
    for w in wins:
        want.setdefault(w["chrom"], []).append((w["start"], w["end"]))
    gw = {}
    with open(GWAS) as fh:
        next(fh)
        for line in fh:
            f = line.split("\t")
            for s, e in want.get(f[0], ()):
                if s < int(f[2]) <= e:
                    gw[(f[0], int(f[2]))] = (float(f[6]), -math.log10(float(f[9])))
                    break
    rows = []
    for w in wins:
        ch, s, e = w["chrom"], w["start"], w["end"]
        geno = {g: obs23a.genotypes(vcf, files[g], ch, s, e) for g in files}
        pos = geno["BAM"][0]
        assert all(np.array_equal(pos, geno[g][0]) for g in geno), "position mismatch between sample sets"
        af = {}
        for g, (_, G) in geno.items():
            called = (G >= 0)
            with np.errstate(invalid="ignore", divide="ignore"):
                af[g] = np.where(called, G, 0).sum(1) / (2.0 * called.sum(1))
        d = af["BAM"] - af["TP1"]
        diff = np.abs(d) >= D_MIN
        bir = np.where(d > 0, af["BIR"], 1 - af["BIR"])[diff]          # Biroundou frequency of the Bambomo-high allele
        hits = [gw[(ch, int(p))] for p in pos[diff] if (ch, int(p)) in gw]
        passing = [lp for a, lp in hits if a > 0.15]
        G = geno["BAM"][1]
        M = (G >= 0).astype(np.int8)
        X = np.where(G >= 0, G, 0).astype(np.int8)
        keep = L.filter_snps(X, M)
        f = frac.get((ch, ((s + e) // 2) // 1_000_000), float("nan"))
        n, t = L.binned_r2(pos[keep], X[keep], M[keep], f)
        r2 = t[LD_BINS].sum() / n[LD_BINS].sum() if n[LD_BINS].sum() else float("nan")
        rows.append(dict(window=w["window"], set=w["set"], window_count=w["count"], n_sites=w["n_sites"],
                         callable_frac=round(f, 4), n_diff=int(diff.sum()), n_diff_in_gwas=len(hits),
                         n_diff_pass_af=len(passing), best_log10p=round(max(passing), 2) if passing else "",
                         bir_mean_freq=round(float(np.nanmean(bir)), 3) if bir.size else "",
                         bir_share_ge025=round(float(np.mean(bir >= 0.25)), 3) if bir.size else "",
                         ld_r2_10_50kb=round(r2, 4), ld_pct_vs_genome=round(100 * float(np.mean(genome <= r2)), 1)))
        print(rows[-1], flush=True)
    with open(out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    main()
