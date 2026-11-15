#!/bin/bash
## 26a_extract_temporal.sh -- genotypes for the temporal analysis (26b, 26c).
## Run from the repository root (reads input_data/ by relative path).
## Groups by collection year from the sample metadata: BAM 2009, BAM 2019 (minus BAM_175, a second-degree
## relative of BAM_167 by KING-robust kinship), APA 2009, APA 2019. Regions: 2 kb every 20 kb on chr1-25, excluding the six
## extended peaks +/- 250 kb and the chr4 inversion. PASS biallelic SNPs; genotypes with FMT/DP < 8 or
## GQ < 20 set missing. Output: one gzipped table per chromosome (CHROM, POS, GT per sample).
##
## Usage (repo root): bash code/11_Drift_vs_Selection/simulation/26a_extract_temporal.sh <vcf> <out_dir> [jobs]
set -euo pipefail
VCF=$1; OUT=$2; JOBS=${3:-6}
mkdir -p "${OUT}/geno"

python3 - "${OUT}" <<'EOF'
import csv, sys
out = sys.argv[1]
drop = {"BAM_175"}
groups = {}
for r in csv.DictReader(open("input_data/01_Terra/data_model/fish_data_2023.txt"), delimiter="\t"):
    pop, sid = r["POP"], r["entity:participant_id"]
    if pop not in ("BAM", "APA") or sid in drop:
        continue
    y = int(r["Date"].split("/")[-1]); y = y + 2000 if y < 100 else y
    groups.setdefault(f"{pop}{y}", []).append(sid)
with open(f"{out}/groups.tsv", "w") as fh:
    for g, ids in sorted(groups.items()):
        for s in ids:
            fh.write(f"{s}\t{g}\n")
        print(g, len(ids))
with open(f"{out}/samples.txt", "w") as fh:
    fh.write("".join(f"{s}\n" for g in sorted(groups) for s in groups[g]))
excl = [("chr4", 40_880_000, 42_770_000)]
for line in open("input_data/10_Introgression/peaks.bed"):
    f = line.split()
    excl.append((f"chr{f[0]}", int(f[1]) - 250_000, int(f[2]) + 250_000))
L = {l.split()[0]: int(l.split()[1]) for l in open(
    "input_data/00_Reference_Genome/jordan-mic4273-mb-hirise-20xx3__06-13-2024__final_assembly.fasta.renamed.fa.fai")}
n = 0
for c in range(1, 26):
    ch = f"chr{c}"
    with open(f"{out}/regions_{ch}.tsv", "w") as fh:
        for s in range(0, L[ch] - 2_000, 20_000):
            e = s + 2_000
            if any(ch == x and s < b and e > a for x, a, b in excl):
                continue
            fh.write(f"{ch}\t{s + 1}\t{e}\n"); n += 1
print("regions", n)
EOF

run_chrom() {
  local ch=$1
  bcftools view -R "${OUT}/regions_${ch}.tsv" -S "${OUT}/samples.txt" -f PASS -v snps -m2 -M2 -Ou "${VCF}" 2>/dev/null \
    | bcftools filter -S . -e 'FMT/DP<8 | FMT/GQ<20' -Ou 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS[\t%GT]\n' | gzip > "${OUT}/geno/${ch}.tsv.gz"
  echo "${ch} $(zcat < "${OUT}/geno/${ch}.tsv.gz" | wc -l)"
}
export -f run_chrom; export VCF OUT
seq 1 25 | sed 's/^/chr/' | xargs -P "${JOBS}" -I{} bash -c 'run_chrom {}'
