#!/bin/bash
## 02_extract_candidates.sh -- candidate ancestry-informative SNPs from the phased VCFs.
## Merges the APA, BAM and BEN phased VCFs (SHAPEIT4, same sites) per chromosome, fills per-panel allele
## frequencies (groups from 01_panels.py: PB/PA primary, FB/FA full, BN Bengue), and keeps every site with
## |AF_BAM - AF_TP| >= 0.5 under any of the three panel pairs (PB-PA, FB-FA, PB-BN). Thinning and pruning per
## panel happen in 03_markers.py. Output per chromosome: cand/<chr>.tsv.gz with CHROM, POS and the phased GT of
## every sample (order in cand/samples.txt). Whole chromosomes are kept; regions are masked later.
##
## Usage (repo root): bash code/11_Drift_vs_Selection/local_ancestry/02_extract_candidates.sh <out_dir> [jobs]
set -euo pipefail
OUT=$1; JOBS=${2:-8}
V=output_data/11_Drift_vs_Selection/sweep_detection
mkdir -p "${OUT}/cand"
bcftools merge "$V/APA.polarized.vcf.gz" "$V/BAM.polarized.vcf.gz" "$V/BEN.polarized.vcf.gz" -r chr25:1-1 -Ou 2>/dev/null \
  | bcftools query -l > "${OUT}/cand/samples.txt"

run_chrom() {
  local ch=$1
  bcftools merge -r "${ch}" "$V/APA.polarized.vcf.gz" "$V/BAM.polarized.vcf.gz" "$V/BEN.polarized.vcf.gz" -Ou 2>/dev/null \
    | bcftools +fill-tags -Ou -- -S "${OUT}/groups_fill.tsv" -t AF 2>/dev/null \
    | bcftools view -i 'abs(INFO/AF_PB-INFO/AF_PA)>=0.5 || abs(INFO/AF_FB-INFO/AF_FA)>=0.5 || abs(INFO/AF_PB-INFO/AF_BN)>=0.5' -Ou \
    | bcftools query -f '%CHROM\t%POS[\t%GT]\n' | gzip > "${OUT}/cand/${ch}.tsv.gz"
  echo "${ch} $(gzip -dc "${OUT}/cand/${ch}.tsv.gz" | wc -l)"
}
export -f run_chrom; export V OUT
seq 1 25 | sed 's/^/chr/' | xargs -P "${JOBS}" -I{} bash -c 'run_chrom {}'
