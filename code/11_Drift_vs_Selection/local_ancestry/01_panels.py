#!/usr/bin/env python3
"""01_panels.py -- sample groups, reference panels and hold-outs for the local-ancestry analysis.

Groups come from the sample metadata (POP, Date, Phenotype), never from the VCF file a sample sits in or its
name prefix: APA_223 is a Bambomo 2019 fish (POP BAM, site 2019-014) that sits in the APA phased VCF.
BAM_175 is dropped (second-degree relative of BAM_167 in a KING-robust kinship check). Clean Apassa TP = Apassa
2009 TP fish with projection score q < 0.05 (26_ancestry_check.tsv, from simulation/26c_ancestry_check.py). Hold-outs are drawn once with seed 20261103:
10 Bambomo 2009, 4 clean Apassa TP, 4 Bengue 2023 TP.

Writes, in <out_dir>:
  samples.tsv      sample, vcf (file it sits in), group, q (26c score or NA), site, roles (comma list)
  groups_fill.tsv  sample -> comma list of panel groups, for bcftools +fill-tags -S
Panel groups: PB/PA primary BAM/APA; FB/FA full BAM/APA; BN Bengue TP reference.

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/01_panels.py <out_dir>
"""
import csv, os, subprocess, sys
import numpy as np

SEED = 20261103
DROP = {"BAM_175"}
META = "input_data/01_Terra/data_model/fish_data_2023.txt"
Q26 = "output_data/11_Drift_vs_Selection/simulation/26_ancestry_check.tsv"
VCFD = "output_data/11_Drift_vs_Selection/sweep_detection"


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    invcf = {}
    for p in ("APA", "BAM", "BEN"):
        names = subprocess.run(["bcftools", "query", "-l", f"{VCFD}/{p}.polarized.vcf.gz"],
                               capture_output=True, text=True, check=True).stdout.split()
        for s in names:
            invcf[s] = p
    q = {}
    for line in open(Q26):
        f = line.rstrip("\n").split("\t")
        if f[0] == "fish":
            q[f[2]] = float(f[3])
    rows = {}
    for r in csv.DictReader(open(META), delimiter="\t"):
        sid, pop = r["entity:participant_id"], r["POP"]
        if pop not in ("APA", "BAM", "BEN") or sid in DROP:
            continue
        y = int(r["Date"].split("/")[-1]); y = y + 2000 if y < 100 else y
        rows[sid] = dict(pop=pop, year=y, pheno=r["Phenotype"], site=r["Site"])
    missing = [s for s in rows if s not in invcf]
    assert not missing, f"in metadata but not in a phased VCF: {missing}"
    extra = [s for s in invcf if s not in rows and s not in DROP]
    assert not extra, f"in a phased VCF but not in metadata groups: {extra}"

    def grp(s):
        r = rows[s]
        return f"{r['pop']}{r['year']}_{r['pheno']}"

    bam09 = sorted(s for s in rows if rows[s]["pop"] == "BAM" and rows[s]["year"] == 2009)
    apa_tp = sorted(s for s in rows if rows[s]["pop"] == "APA" and rows[s]["year"] == 2009 and rows[s]["pheno"] == "TP")
    clean = [s for s in apa_tp if q[s] < 0.05]
    ben = sorted(s for s in rows if rows[s]["pop"] == "BEN" and rows[s]["pheno"] == "TP")
    rng = np.random.default_rng(SEED)
    ho_bam = set(rng.choice(bam09, 10, replace=False))
    ho_apa = set(rng.choice(clean, 4, replace=False))
    ho_ben = set(rng.choice(ben, 4, replace=False))

    roles = {s: [] for s in rows}
    for s in bam09:
        roles[s].append("FB")
        roles[s].append("holdout_BAM" if s in ho_bam else "PB")
    for s in clean:
        roles[s].append("FA")
        roles[s].append("holdout_APA" if s in ho_apa else "PA")
    for s in ben:
        roles[s].append("holdout_BEN" if s in ho_ben else "BN")
    with open(f"{out}/samples.tsv", "w") as fh:
        fh.write("sample\tvcf\tgroup\tq26c\tsite\troles\n")
        for s in sorted(rows, key=lambda x: (grp(x), x)):
            fh.write(f"{s}\t{invcf[s]}\t{grp(s)}\t{q.get(s, float('nan')):.3f}\t{rows[s]['site']}\t"
                     f"{','.join(roles[s]) or '-'}\n")
    with open(f"{out}/groups_fill.tsv", "w") as fh:
        for s in sorted(rows):
            g = [x for x in roles[s] if x in ("PB", "PA", "FB", "FA", "BN")]
            if g:
                fh.write(f"{s}\t{','.join(g)}\n")
    for g in ("PB", "PA", "FB", "FA", "BN", "holdout_BAM", "holdout_APA", "holdout_BEN"):
        mem = sorted(s for s in rows if g in roles[s])
        print(g, len(mem), " ".join(mem) if len(mem) <= 13 else "")
    from collections import Counter
    print(Counter(grp(s) for s in rows))


if __name__ == "__main__":
    main()
