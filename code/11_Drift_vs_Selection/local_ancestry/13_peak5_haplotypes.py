#!/usr/bin/env python3
"""13_peak5_haplotypes.py -- what are the Bambomo-like Peak 5 copies in Apassa? (post hoc)

12_peak_runs.py finds Bambomo-like Peak 5 copies in Apassa (APA_6681 and APA_198 on both haplotypes), although
no Apassa fish carries the BP allele at the Peak 5 index SNP. Question: do those copies resemble
Bambomo haplotypes that carry the TP index allele (as BAM_6597, the Bambomo TP fish, does on both copies), Bambomo
BP-allele haplotypes, or other Apassa haplotypes? Descriptive similarity only (similar is not the same).

Data: the per-population phased VCFs used by 02_extract_candidates.sh (sweep_detection/APA and BAM .polarized.vcf.gz, SHAPEIT4; same haplotype order as
hap_primary), extended Peak 5 chr17:11,494,062-11,781,225. Sites: biallelic, minor-allele count >= 2 over all APA
and BAM haplotypes. Distance: fraction of sites that differ between two haplotypes. Index allele at
chr17:11,692,887, oriented so BP = the majority allele among BAM 2009 BP haplotypes.
Reference sets: BAM haplotypes with the BP index allele; BAM haplotypes with the TP index allele (BAM_6597 listed
separately); other Apassa haplotypes (not in results/peak_runs.tsv for Peak 5; this includes the Apassa panel
fish, whose posteriors were not computed, so it may hold a few Bambomo-like copies). Focal haplotypes: every
Bambomo-like Apassa Peak 5 copy, plus both BAM_6597 haplotypes.
For each: distance to each set (minimum and median, other fish only) and the five nearest haplotypes from other
fish.

Writes <out_dir>/results/peak5_haplotypes.tsv and results/peak5_haplotypes_nearest.tsv.

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/13_peak5_haplotypes.py <out_dir>
"""
import csv, subprocess, sys
import numpy as np

V = "output_data/11_Drift_vs_Selection/sweep_detection"
REG = "chr17:11494062-11781225"
IDX = 11692887


def main():
    O = sys.argv[1]
    info = {}
    for line in open(f"{O}/samples.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "sample":
            info[f[0]] = f[2]
    cmd = f"bcftools merge -r {REG} {V}/APA.polarized.vcf.gz {V}/BAM.polarized.vcf.gz -Ou 2>/dev/null"
    names = subprocess.run(cmd + " | bcftools query -l", shell=True, capture_output=True, text=True).stdout.split()
    out = subprocess.run(cmd + " | bcftools view -m2 -M2 -v snps -Ou | bcftools query -f '%POS[\\t%GT]\\n'",
                         shell=True, capture_output=True, text=True).stdout.splitlines()
    pos, rows = [], []
    for line in out:
        f = line.split("\t")
        g = []
        for x in f[1:]:
            if "|" not in x or "." in x:
                g = None
                break
            a, b = x.split("|")
            g += [int(a), int(b)]
        if g is not None:
            pos.append(int(f[0])); rows.append(g)
    pos = np.array(pos); H = np.array(rows, np.int8)
    hap_fish = np.repeat(names, 2); hap_k = np.tile([0, 1], len(names))
    keep = [i for i, s in enumerate(hap_fish) if s in info and (info[s].startswith("APA") or info[s].startswith("BAM"))]
    H, hap_fish, hap_k = H[:, keep], hap_fish[keep], hap_k[keep]
    ac = H.sum(1)
    ok = (ac >= 2) & (ac <= H.shape[1] - 2)
    j = np.flatnonzero(pos == IDX)
    assert j.size == 1, "index SNP missing"
    idx_allele = H[j[0]].copy()
    H = H[ok]
    grp = np.array([info[s] for s in hap_fish])
    bp_maj = round(idx_allele[grp == "BAM2009_BP"].mean())
    is_bp = idx_allele == bp_maj
    car = set()
    for r in csv.DictReader(open(f"{O}/results/peak_runs.tsv"), delimiter="\t"):
        if r["peak"] == "5":
            car.add((r["sample"], int(r["hap"])))
    focal = [i for i in range(len(hap_fish)) if (hap_fish[i], hap_k[i]) in car or hap_fish[i] == "BAM_6597"]
    bam = np.array([g.startswith("BAM") for g in grp])
    is6597 = hap_fish == "BAM_6597"
    sets = {"BAM_BP_index_allele": bam & is_bp,
            "BAM_TP_index_allele_excl_6597": bam & ~is_bp & ~is6597,
            "BAM_6597": is6597,
            "APA_other": ~bam & np.array([(s, k) not in car for s, k in zip(hap_fish, hap_k)])}
    D = (H[:, focal, None] != H[:, None, :]).mean(0)          # focal x all
    res = [f"# {ok.sum()} sites (MAC >= 2) in {REG}; BP index allele code {bp_maj} (polarized VCF); "
           + "; ".join(f"{k} n={v.sum()}" for k, v in sets.items()),
           "sample\thap\tgroup\tindex_allele\t" + "\t".join(f"{k}_min\t{k}_median" for k in sets)]
    near = ["sample\thap\trank\tother\tother_hap\tother_group\tother_index_allele\tdistance"]
    for a, i in enumerate(focal):
        d = D[a]
        other = hap_fish != hap_fish[i]
        line = f"{hap_fish[i]}\t{hap_k[i]}\t{grp[i]}\t{'BP' if is_bp[i] else 'TP'}"
        for k, m in sets.items():
            v = d[m & other]
            line += f"\t{v.min():.4f}\t{np.median(v):.4f}" if v.size else "\tNA\tNA"
        res.append(line)
        cand = np.flatnonzero(other)
        for r, b in enumerate(cand[np.argsort(d[cand], kind="stable")][:5], 1):
            near.append(f"{hap_fish[i]}\t{hap_k[i]}\t{r}\t{hap_fish[b]}\t{hap_k[b]}\t{grp[b]}\t"
                        f"{'BP' if is_bp[b] else 'TP'}\t{d[b]:.4f}")
    open(f"{O}/results/peak5_haplotypes.tsv", "w").write("\n".join(res) + "\n")
    open(f"{O}/results/peak5_haplotypes_nearest.tsv", "w").write("\n".join(near) + "\n")
    print("\n".join(res)); print("\n".join(near))


if __name__ == "__main__":
    main()
