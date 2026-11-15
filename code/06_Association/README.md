# 06_Association: GEMMA genome-wide association of BP vs TP EODs

Produces: Fig. 3A (`gwas_report.Rmd`), Table S1 (`make_table_s1.R`), the permutation thresholds, peak definitions and
Wald vs likelihood-ratio comparison of the ESM Association Analysis, and the Bambomo mixed-vs-BP scan of the
Supplementary Note (prepared here, run in `14_Wobble_Analysis`). The peak files feed 08, 09, 10, 11, 14 and 15.

## Inputs
- `input_data/06_Association/<COMPARISON>.keep.txt` / `.pheno.txt` — sample list and phenotype code per comparison
  (shipped; see "Keep and phenotype files" below):
  - `PKINGS_ALL_WOB_EXCLUDED` — main GWAS: all 240 *P. kingsleyae*, TP = 1, BP = 2, mixed = -9 (198 analysed).
    Also read by `03_VCF_QC_Filter/ref_bias_summary.R`.
  - `BAM_ONLY_BP0_TPWOB1` — Bambomo only (67 fish): BP = 0, TP and mixed = 1. Its prepared genotypes and kinship are
    reused by the mixed-vs-BP scans and power simulation in `14_Wobble_Analysis` (07–11).
  - `PKINGS_TP1_BP2_WOB9` — read only by `13_Structure/prep_admixture.sb`. Same 240 fish and codes as
    `PKINGS_ALL_WOB_EXCLUDED` except that BAM_6494, BAM_6605, BIK_6908 and BIR_344 are coded -9.
- `input_data/06_Association/PKINGS_ALL_WOB_EXCLUDED_GENES_IN_PEAKS_manual_annos_add.txt` — hand-curated gene
  names used by `make_table_s1.R` after `input_data/08_Peak_Analysis/pretty_gene_names.csv`.
- `output_data/03_QC/${snps_only_vcf}` — filtered biallelic SNPs (03_VCF_QC_Filter).
- `input_data/00_Reference_Genome/paramormyrops_ncbi_0.4.renamed.gff` — RefSeq annotation GCF_048594095.1-RS_2025_03
  with chromosomes renamed (see the reference-genome section of the top-level README); read by `make_table_s1.R`.
- Containers: `${gwas_tools_image}` (bcftools, qctool, GEMMA, PLINK 1.9, R). Permutations also need
  [gemma-wrapper](https://github.com/genetics-statistics/gemma-wrapper) 0.99.1, which `04_gemma_permutation_local.sh`
  expects at `~/gemma-wrapper-0.99.1/` and runs in `--slurm` mode (edit the path for your installation).

## Keep and phenotype files
The shipped files are the ones used for the paper; use them to reproduce its numbers. Row order sets the sample order of
the GEMMA genotype and phenotype files and therefore the permuted phenotypes. Regenerating them from the sample sheet
gives the same samples and codes, but in sample-sheet order.

They were written by `generate_phenotype_and_keepfiles.py`, which asks for a code for each phenotype present
(`TP`, `BP`, `Wobble` = mixed) and then for the output name:

```bash
source config.env
Rscript code/00_common/make_sample_sheets.R          # writes input_data/01_Terra/data_model/fish_data_*.txt from Data S1
kfiledir=${root}/input_data/06_Association
gen() { singularity exec --bind $root:/project_root --bind $kfiledir:/out_dir ${gwas_tools_image} \
  python3 /project_root/code/06_Association/generate_phenotype_and_keepfiles.py \
  /project_root/input_data/01_Terra/data_model/fish_data_2024a.txt --output_directory /out_dir "$@"; }

gen --include_species "Paramormyrops kingsleyae"   # TP 1, BP 2, Wobble -9; name PKINGS_ALL_WOB_EXCLUDED
gen --include_pop BAM                              # TP 1, BP 0, Wobble 1;  name BAM_ONLY_BP0_TPWOB1
```

`PKINGS_TP1_BP2_WOB9` was written the same way as `PKINGS_ALL_WOB_EXCLUDED`, from a phenotype table that coded the
four fish listed above as missing.

## Steps
Run from the repository root after `source config.env`. `<C>` is a comparison name. The `.sb` scripts are SLURM jobs;
edit their `#SBATCH` lines (account, partition, resources) for your cluster.

1. `bash code/06_Association/01_run_prep_data.sh <C>` — subset the SNP VCF to the keep file and drop SNPs with
   MAF ≤ 0.05, > 10% missing data or monomorphic (`prep_data.sb`); writes `<C>_merge_ID.vcf.gz` and the SNP map.
2. `bash code/06_Association/02_run_rel_mat.sh <C>` — thin to one random SNP per 20 kb (seed 42), write the
   BIMBAM dosages and the centred relatedness matrix (`gemma -gk 1`; `relate_matrix.sb`).
3. `bash code/06_Association/03_run_gemma.sh <C>` — GEMMA `-lmm 2` (likelihood-ratio test) on all SNPs
   (`gemma.sb`); writes `<C>.lmm2.assoc.txt`.
4. `. code/06_Association/04_gemma_permutation_local.sh <C>` — 100 phenotype permutations with gemma-wrapper.
5. `. code/06_Association/05_collect_min_values.sh <C>` — minimum p of each permutation into
   `<C>_permution/min.txt` (95th percentile = significant, 67th = suggestive).
6. `bash code/06_Association/06_Render_GWAS_report.sh <C> 0.15` — render `gwas_report.Rmd` (thresholds, alternate
   allele frequency > 0.15 filter, peak calling, Manhattan plot). Optional third argument `lrt` (default) or `wald`.
7. `bash code/06_Association/07_get_ld_snps.sh <C>` — PLINK r² between each peak SNP and SNPs within 10,000 SNPs or
   10 Mb (`<C>_index_snps.ld.ld`, read by 08).
8. `Rscript code/06_Association/make_table_s1.R` — Table S1 (after 08 `peak_analysis_new.Rmd` has written
   `PKINGS_ALL_WOB_EXCLUDED_GENES_IN_PEAKS.txt`).

For the main GWAS run steps 1–7 with `<C> = PKINGS_ALL_WOB_EXCLUDED`. For the Bambomo scan run steps 1–2 with
`BAM_ONLY_BP0_TPWOB1` on the phased call set (the `prep_data.sb` settings of the Wald comparison below) and continue
in `14_Wobble_Analysis`.

### Core peaks and extended windows
`gwas_report.Rmd` merges af-filtered SNPs with −log10 p > 3 within 75 kb (`peak_list_permutation(pk.lmm.filt, "log_p",
3, 4)`) and keeps regions with > 5 SNPs above the suggestive threshold. Rendered as committed, it gives the extended
windows: that run's `PEAKS.bed` is `input_data/10_Introgression/peaks.bed` (shipped), and its 645 peak SNPs
(`_SNPS_IN_PEAKS.txt`, `_SNPS_IN_PEAKS_TOP250.txt`) are the SNP lists used by step 7 and by
`08_Peak_Analysis/01_get_SNP_HITS.sh`. The core peaks of Table S1 (`PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed`, read by
`make_table_s1.R`, 08 and 09) come from rendering the same notebook with that merge threshold set to 5,
which overwrites `PEAKS.bed`, `_SNPS_IN_PEAKS*.txt` and `_index_snps.txt` in the same directory. Both definitions give
the same six peaks and index SNPs.

### Wald vs likelihood-ratio comparison
The ESM compares a Wald run and a matched likelihood-ratio run on the phased call set (07_Phasing). The Wald run is
steps 1–6 for `PKINGS_ALL_WOB_EXCLUDED` with `indir=${root}/output_data/07_Phasing/` and
`snps_only_vcf=${phased_vcf}` in `prep_data.sb`, `-lmm 1` in place of `-lmm 2` in `gemma.sb` and
`04_gemma_permutation_local.sh` (output `PKINGS_ALL_WOB_EXCLUDED.lmm1.assoc.txt`), and the report rendered with
`wald`. It writes to the same directory as the main run, so run it in a separate copy of the repository and first copy
the main run's `PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed` and `_permution/min.txt` to
`output_data/06_Association/WALD_VS_LRT/published_PEAKS.bed` and `published_min.txt`. The matched run applies
GEMMA `-lmm 2` to the same phased genotypes, kinship and 100 permuted phenotypes, in
`output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED_LRT/` (files `PKINGS_ALL_WOB_EXCLUDED_LRT.*`, report rendered with
`lrt`). Then `sbatch code/06_Association/compare_wald_lrt.sb` runs `compare_wald_lrt.R`.

## Outputs
- `output_data/06_Association/<C>/<C>.lmm2.assoc.txt` — GEMMA results.
- `output_data/06_Association/<C>/<C>_permution/min.txt` — permutation minima.
- `output_data/06_Association/<C>/<C>_PEAKS.bed`, `<C>_SNPS_IN_PEAKS.txt`, `<C>_SNPS_IN_PEAKS_TOP250.txt`,
  `<C>_index_snps.txt`, `<C>.gwas.svg` (Fig. 3A), `<C>_GWAS_Report.report.nb.html`.
- `output_data/06_Association/<C>/<C>_index_snps.ld.ld` — LD with the peak SNPs.
- `output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/Table_S1_core_peaks.tsv`, `Table_S1_core_peak_genes.tsv`.
- `output_data/06_Association/WALD_VS_LRT/` — `summary.tsv`, `published_peaks_check.tsv`, `snp_agreement.tsv`,
  `perm_min_paired.tsv`, `new_peaks_not_in_published.tsv` and two PNG plots.
