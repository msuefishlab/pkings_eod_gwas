# 06_peak6_curv3: Peak 6 ABBA-BABA with the *curvifrons* group (CURV3)

Produces: the "Targeted test at Peak 6" of the ESM Introgression Methods and its rows in Supplementary Data S7.
Tests whether the Peak 6 BP2 haplotype is shared with *P. curvifrons*, *P. ntotom* and *P.* sp. 'SN2' (CURV3) beyond
the genome-wide level, with the other groups, peak interval and Dsuite version of the main analysis.

## Inputs
- `input_data/10_Introgression/SETS_curv3.txt` — `SETS.txt` with the 5 PCUR, 4 PNTO and 4 PSN2 fish as `CURV3`, and
  the 3 PLON fish and PSZA_455 as `xxx` (Dsuite drops `xxx`). All other assignments unchanged.
- `input_data/10_Introgression/pops_curv3.tree` — `((((((BP1,TP1),((BP2,TP2),BP3)),TP_OUT),CURV3),Outgroup));`
- `input_data/10_Introgression/test_trios_curv3.txt` — the six trios for Dinvestigate (P1 = lineage tested,
  P2 = its reference, P3 = CURV3).
- `input_data/10_Introgression/test_trios_published_check.txt` — BP2 TP_OUT BP_ANC, to check the chr24 windows
  against the main analysis.
- `input_data/10_Introgression/SETS.txt`, `pops.tree`, `peaks.bed` — the main analysis's groups, tree and extended
  peak windows.
- `output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz` — the
  biallelic SNP shards plus the indel shards (`code/03_VCF_QC_Filter/05_combine_filtered.sh`); `bcftools view -v snps
  -m2 -M2` gives back the SNP-only VCF of the main analysis. Pass another path as the first argument of 01–03.
- Main-analysis outputs in `output_data/10_Introgression/` (parent README): `whole_genome_combined_tree.txt`,
  `peaks_extended/peak_6_tree.txt`, `peak_fdm_summary.tsv` and `BP2_TP_OUT_BP_ANC_localFstats_chr24_50_25.txt`.
- `input_data/10_Introgression/SupplementalDataS7_peak_ABBA_BABA.tsv` — Supplementary Data S7, downloaded from the paper.
- Docker, bcftools, tabix and R (tidyverse).

## Dsuite version
Dsuite 0.5 r58, as in the main analysis: commit a547f99, whose only change after the r58 commit (8b0bb96) is an
`#include`. The project Dockerfile (`images/DsuiteDocker`) builds it (image `dsuite_local`). On aarch64 (e.g. Apple
silicon) that binary hangs: Dsuite reads options with `for (char c; (c = getopt_long(...)) != -1;)`, and `char` is
unsigned on ARM, so the loop never ends. `Dockerfile.dsuite_arm` recompiles the same commit with `-fsigned-char`
(x86 semantics; no source change; image `dsuite_arm`). The scripts call `dsuite_arm`.

## Steps
Run from `code/10_Introgression/06_peak6_curv3/`.
1. `bash 00_build_dsuite.sh` — build the Docker images and record the version (`dsuite_version.txt`).
2. `bash 01_validate_peak6.sh` — rerun the main trio set on Peak 6 and compare with Data S7 (`validate_vs_dataS7.R`).
   CURV3 results are read only if this reproduces Data S7 (BP2–BP_ANC D = 0.197, 445 informative sites).
3. `bash 02_make_chr_vcfs.sh "" 12` — chr1–25 biallelic SNPs, the 209 SETS samples, GT only (`_scratch/`).
4. `bash 03_run_dtrios.sh "" 8` — Peak 6 with the CURV3 sets; genome-wide CURV3 and main sets per chromosome;
   DtriosCombine.
5. `bash 04_run_dinvestigate.sh 7` — Dinvestigate `-w 50,25` per chromosome for the six CURV3 trios, and the main
   BP2–TP_OUT–BP_ANC trio on chr24.
6. `Rscript validate_vs_dataS7.R` — again, now including the genome-wide rerun of the main sets.
7. `Rscript 05_summarise.R` — summary table, f_dM block percentiles and verdict.
8. `bash 06_cleanup.sh` — delete `_scratch/` and gzip the Dinvestigate windows and logs.

Steps 4 and 5 can run at the same time. The Dsuite calls are those of `code/10_Introgression/run_dsuite*.sb`, with
Docker in place of Singularity:

```bash
bcftools view -r chr24:8331219-8730509 -v snps -m2 -M2 $VCF -Oz -o peak_6.vcf.gz
Dsuite Dtrios -o peak_curv3/peak_6 -t pops_curv3.tree peak_6.vcf.gz SETS_curv3.txt
bcftools view -r chrN -v snps -m2 -M2 -S sets_samples.txt -Ou $VCF | bcftools annotate -x INFO,^FORMAT/GT -Oz -o chrN.snps.gt.vcf.gz
Dsuite Dtrios -o genome_curv3/chrN -t pops_curv3.tree chrN.snps.gt.vcf.gz SETS_curv3.txt
Dsuite DtriosCombine -o genome_curv3/whole_genome -t pops_curv3.tree genome_curv3/chr1 ... genome_curv3/chr25
Dsuite Dinvestigate -w 50,25 -n chrN chrN.snps.gt.vcf.gz SETS_curv3.txt test_trios_curv3.txt
```

Dropping samples outside the SETS file, INFO and FORMAT fields other than GT keeps every record, so the Dsuite input
(and the jackknife blocks, which are set by the number of variant lines) is unchanged. The f4-ratio differs slightly
between runs (≤ 0.007) because Dsuite splits P3 at random on every run; D, Z and the site-pattern counts are
deterministic.

## Summary statistics
`05_summarise.R` reuses the main analysis's code: `orient_trios()` from `DsuiteAnalysis.Rmd` (D signed so + = sharing
between the named lineage and CURV3; Z_vs_genome = (D_peak − D_genome)/SE_peak, SE_peak = |D/Z|), and
`summarize_rep_side_globalbg()` (`dinvestigate_functions.R`) plus `build_peak_summary()`, read from
`DsuiteInvestigateAnalysis.Rmd` at run time (f_dM × −1, peak median over windows overlapping the extended interval,
background = windows outside all six peaks, block null = medians of runs of the same number of consecutive background
windows). q is BH over the six tests. `signal` is the full rule; `D_signal` is its D part (the `signal` column of
Data S7). A TP control counts as "showing the same excess" if it meets the full rule.

## Outputs
`output_data/10_Introgression/peak6_curv3/`:
- `peak6_curv3_dtrios_summary.tsv`, `rule_verdict.txt` — results.
- `validation_vs_dataS7.tsv`, `validation_peak6_all_trios_vs_published.tsv`,
  `validation_genome_published_sets_vs_published.tsv`, `check_fdm_code_vs_published_peak_fdm_summary.tsv`,
  `check_dinvestigate_chr24_vs_published.tsv` — reproduction checks.
- `peak_curv3/`, `genome_curv3/`, `genome_published_sets/`, `validation_published_sets/`, `dinvestigate_curv3/`,
  `dinvestigate_published_check/`, `logs/`, `dsuite_version.txt`.
