#!/usr/bin/env python3
"""14_peak5_diag_snps.py -- Peak 5 beyond the index SNP: is any part of the Apassa copies a Bambomo BP block? (post hoc)

Checks the Peak 5 reading of 12/13 beyond the index SNP alone. This uses the haplotype-audit SNP set and reference
pools (code/14_Wobble_Analysis/17_peak_hap_match_audit.R, section 2):
the Fig. S4 Peak 5 SNPs (TOP250 of the extended window, chr17), genotypes from the 07_Phasing phased VCF, BP pool =
Bambomo BP fish, TP pool = APA/BEN/DOV/DOG TP fish, consensus = majority allele per pool, diagnostic = consensus
differs. For every haplotype of the Apassa, Bengue and Bambomo fish:
  bp_string          BP (1) or TP (0) allele at each diagnostic SNP in position order ('|' marks the index SNP)
  n_bp / n_diag      BP alleles carried
  longest_bp_run     longest run of consecutive diagnostic SNPs carrying the BP allele, its span and whether it
                     contains the index SNP
  n_bp_left/right    BP alleles at diagnostic SNPs left of the index SNP / at and right of it
  nn_mismatch        mismatches over all window SNPs to the nearest Bambomo BP-class copy (other fish only)
  longest_ident_run  longest run of consecutive window SNPs identical to some Bambomo BP-class copy (other fish
                     only), its span and whether it contains the index SNP
Bambomo BP-class copy: a haplotype of a Bambomo BP fish that is nearer the BP than the TP consensus (n_bp > half
the diagnostic SNPs), the 01/17 two-class rule, as the audit's home copies.
Calibration rows are every haplotype, so Bambomo BP copies (against other BP fish) and Apassa copies that are not
Bambomo-like show what the statistics look like for "same" and for "unrelated".

Per fish (Apassa, Bengue, Bambomo): copies carrying the left BP block (>= 90% of the left diagnostic SNPs BP) and
the right BP block (>= 90% of the at/right SNPs BP); 90% is the audit's full/partial tier boundary, applied per half
after the halves were seen (post hoc). Counts by population x phenotype x (left, right) genotype, and an exact
Fisher test of left-block homozygotes, Apassa wobble against Apassa TP.

Writes <out_dir>/results/peak5_diag_snps.tsv, results/peak5_diag_snps_summary.tsv and
results/peak5_half_genotypes.tsv.

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/14_peak5_diag_snps.py <out_dir>
"""
import csv, subprocess, sys
import numpy as np

S4 = "output_data/08_Peak_Analysis/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS_TOP250.vcf.gz"
PH = ("output_data/07_Phasing/"
      "MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_SNPS_ONLY.phased.renamed.vcf.gz")
META = "input_data/01_Terra/data_model/fish_data_2023.txt"
CHR, IDX = "chr17", 11692887


def runs(v):
    """Longest run of True in v: (length, start index, end index)."""
    best, cur, s, bs = 0, 0, 0, 0
    for i, x in enumerate(v):
        if x:
            if cur == 0:
                s = i
            cur += 1
            if cur > best:
                best, bs = cur, s
        else:
            cur = 0
    return best, bs, bs + best - 1


def main():
    O = sys.argv[1]
    meta = {r["entity:participant_id"]: r for r in csv.DictReader(open(META), delimiter="\t")}
    sites = subprocess.run(f"bcftools query -f '%CHROM\\t%POS\\n' {S4}", shell=True,   # file has no index
                           capture_output=True, text=True).stdout
    sites = "".join(l + "\n" for l in sites.splitlines() if l.split("\t")[0] == CHR)
    tf = f"{O}/results/.peak5_sites.tsv"
    open(tf, "w").write(sites)
    names = subprocess.run(f"bcftools query -l {PH}", shell=True, capture_output=True, text=True).stdout.split()
    out = subprocess.run(f"bcftools view -T {tf} -R {tf} {PH} 2>/dev/null | bcftools query -f '%POS[\\t%GT]\\n'",
                         shell=True, capture_output=True, text=True).stdout.splitlines()
    pos = np.array([int(l.split("\t")[0]) for l in out])
    assert len(pos) == len(sites.splitlines()), "a Fig. S4 SNP is missing from the phased VCF"
    G = [l.split("\t")[1:] for l in out]
    H = np.full((len(pos), 2 * len(names)), -1, np.int8)
    for i, g in enumerate(G):
        for j, x in enumerate(g):
            if "|" in x and "." not in x:
                a, b = x.split("|"); H[i, 2 * j] = int(a); H[i, 2 * j + 1] = int(b)
    hf = np.repeat(names, 2); hk = np.tile(["A", "B"], len(names))
    pop = np.array([meta[s]["POP"] if s in meta else "" for s in hf])
    phe = np.array([meta[s]["Phenotype"] if s in meta else "" for s in hf])
    bpref = (pop == "BAM") & (phe == "BP")
    tpref = np.isin(pop, ["APA", "BEN", "DOV", "DOG"]) & (phe == "TP")

    def freq(m):
        X = H[:, m].astype(float); X[X < 0] = np.nan
        return np.nanmean(X, 1)
    bf, tfq = freq(bpref), freq(tpref)
    bpc, tpc = (bf >= 0.5).astype(int), (tfq >= 0.5).astype(int)
    diag = np.flatnonzero(bpc != tpc)
    ii = np.flatnonzero(pos == IDX)[0]
    idx_in_diag = ii in diag
    show = np.isin(pop, ["APA", "BEN", "BAM"])
    rows = ["sample\thap\tPOP\tPhenotype\tn_diag\tn_bp\tn_bp_left\tn_bp_right\tlongest_bp_run_snps\tlongest_bp_run_span\t"
            "bp_run_has_index\tnn_mismatch\tnearest_BP_hap\tlongest_ident_run_snps\tlongest_ident_run_span\t"
            "ident_run_has_index\tindex_allele\tbp_string"]
    res = []
    nbp_all = (H[diag] == bpc[diag][:, None]).sum(0)
    bpcols = np.flatnonzero(bpref & (nbp_all > len(diag) / 2))
    left = pos[diag] < IDX
    for c in np.flatnonzero(show):
        h = H[:, c]
        isbp = h[diag] == bpc[diag]
        n, a, b = runs(isbp)
        bpspan = f"{pos[diag][a]}-{pos[diag][b]}" if n else "-"
        bphas = bool(n) and pos[diag][a] <= IDX <= pos[diag][b]
        others = [k for k in bpcols if hf[k] != hf[c]]
        mm = np.array([((h != H[:, k]) & (h >= 0) & (H[:, k] >= 0)).sum() for k in others])
        knn = others[int(np.argmin(mm))]
        best = (0, 0, 0)
        for k in others:
            r = runs((h == H[:, k]) | (h < 0) | (H[:, k] < 0))
            if r[0] > best[0]:
                best = r
        ispan = f"{pos[best[1]]}-{pos[best[2]]}"
        ihas = pos[best[1]] <= IDX <= pos[best[2]]
        s = "".join(str(int(x)) for x in isbp)
        k = int(np.searchsorted(pos[diag], IDX))
        s = s[:k] + "|" + s[k:]
        idx_allele = "BP" if h[ii] == (1 if bf[ii] > tfq[ii] else 0) else "TP"
        nl, nr = int(isbp[left].sum()), int(isbp[~left].sum())
        rec = dict(sample=hf[c], hap=hk[c], POP=pop[c], phe=phe[c], n_bp=int(isbp.sum()), run=n, nn=int(mm.min()),
                   irun=best[0], ihas=ihas, idx=idx_allele)
        rec["L"] = nl >= 0.9 * left.sum(); rec["R"] = nr >= 0.9 * (~left).sum()
        res.append(rec)
        rows.append(f"{hf[c]}\t{hk[c]}\t{pop[c]}\t{phe[c]}\t{len(diag)}\t{rec['n_bp']}\t{nl}\t{nr}\t{n}\t{bpspan}\t{int(bphas)}\t"
                    f"{rec['nn']}\t{hf[knn]}|{hk[knn]}\t{best[0]}\t{ispan}\t{int(ihas)}\t{idx_allele}\t{s}")
    open(f"{O}/results/peak5_diag_snps.tsv", "w").write("\n".join(rows) + "\n")

    car = set()
    for r in csv.DictReader(open(f"{O}/results/peak_runs.tsv"), delimiter="\t"):
        if r["peak"] == "5":
            car.add(r["sample"])
    groups = {"BAM BP fish, BP index allele": lambda r: r["POP"] == "BAM" and r["phe"] == "BP" and r["idx"] == "BP",
              "BAM, TP index allele": lambda r: r["POP"] == "BAM" and r["idx"] == "TP",
              "APA, Bambomo-like at Peak 5 in 12_peak_runs (fish level)": lambda r: r["POP"] == "APA" and r["sample"] in car,
              "APA, other fish": lambda r: r["POP"] == "APA" and r["sample"] not in car,
              "BEN": lambda r: r["POP"] == "BEN"}
    sm = [f"# {len(pos)} window SNPs ({pos.min()}-{pos.max()}), {len(diag)} diagnostic; index SNP chr17:{IDX} "
          f"diagnostic: {idx_in_diag}; diagnostic SNPs left of index {left.sum()} ({pos[diag][left].min()}-"
          f"{pos[diag][left].max()}), at/right {(~left).sum()} ({pos[diag][~left].min()}-{pos[diag][~left].max()}); "
          f"Bambomo BP-class reference copies {len(bpcols)}",
          "group\tn_haps\tn_bp_median\tn_bp_range\tlongest_bp_run_median\tnn_mismatch_median\tnn_mismatch_range\t"
          "longest_ident_run_median\tlongest_ident_run_range\tident_run_has_index"]
    for g, f in groups.items():
        R = [r for r in res if f(r)]
        if not R:
            continue
        v = lambda k: np.array([r[k] for r in R])
        sm.append(f"{g}\t{len(R)}\t{np.median(v('n_bp')):g}\t{v('n_bp').min()}-{v('n_bp').max()}\t"
                  f"{np.median(v('run')):g}\t{np.median(v('nn')):g}\t{v('nn').min()}-{v('nn').max()}\t"
                  f"{np.median(v('irun')):g}\t{v('irun').min()}-{v('irun').max()}\t{int(v('ihas').sum())}")
    from collections import Counter, defaultdict
    from scipy.stats import fisher_exact
    fishg = defaultdict(lambda: [0, 0])
    fp = {}
    for r in res:
        fishg[r["sample"]][0] += int(r["L"]); fishg[r["sample"]][1] += int(r["R"]); fp[r["sample"]] = (r["POP"], r["phe"])
    tab = Counter((fp[f][0], fp[f][1], L, R) for f, (L, R) in fishg.items())
    ht = ["POP\tPhenotype\tleft_block_copies\tright_block_copies\tn_fish\tfish"]
    for k in sorted(tab):
        fs = sorted(f for f in fishg if (fp[f][0], fp[f][1], *fishg[f]) == k)
        ht.append(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{tab[k]}\t{','.join(fs) if len(fs) <= 12 else '-'}")
    aw = [f for f in fishg if fp[f] == ("APA", "Wobble")]; at = [f for f in fishg if fp[f] == ("APA", "TP")]
    hw = sum(fishg[f][0] == 2 for f in aw); htp = sum(fishg[f][0] == 2 for f in at)
    p = fisher_exact([[hw, len(aw) - hw], [htp, len(at) - htp]])[1]
    ht.append(f"# Apassa left-block homozygotes: wobble {hw}/{len(aw)}, TP {htp}/{len(at)}; Fisher exact two-sided "
              f"p = {p:.3f} (post hoc; wobbles were picked at 2019-017)")
    open(f"{O}/results/peak5_half_genotypes.tsv", "w").write("\n".join(ht) + "\n")
    print("\n".join(ht))
    open(f"{O}/results/peak5_diag_snps_summary.tsv", "w").write("\n".join(sm) + "\n")
    print("\n".join(sm))


if __name__ == "__main__":
    main()
