# 08_Peak_Analysis: regional plots, peak genes and variance explained at the GWAS peaks

Produces: Fig. 3B (`peak_analysis_new.Rmd`) and Fig. 3C (`PVE_new.Rmd`), plus the list of genes in or near each peak
(`PKINGS_ALL_WOB_EXCLUDED_GENES_IN_PEAKS.txt`) read by `06_Association/make_table_s1.R` and 09_RNASeq.
`peak_analysis_functions.R` is also sourced by `05_PopGen/peak_popgen*.Rmd`.

## Inputs
- `output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/` — GEMMA results, `min.txt`, core `_PEAKS.bed`,
  `_SNPS_IN_PEAKS*.txt` and `_index_snps.ld.ld` (06_Association; see its README for the core and extended peak runs).
- `output_data/03_QC/${snps_only_vcf}` — unphased SNP VCF (03_VCF_QC_Filter).
- `input_data/04_Phylogenomics/pkings_iqtree_genes_astral_gcf_scf_tree.cf.tree.ultra.time.nex` — species tree (time-calibrated; shipped because it was written interactively from the 04 Paramormyrops notebook)
  (04_Phylogenomics), drawn beside the genotype matrix by `phylogenoplot_multi.R`.
- `input_data/00_Reference_Genome/paramormyrops_ncbi_0.4.renamed.gff` — RefSeq annotation GCF_048594095.1-RS_2025_03
  with chromosomes renamed (top-level README).
- `input_data/01_Terra/data_model/fish_data_2024a.txt` — sample sheet (`Rscript code/00_common/make_sample_sheets.R`).
- `input_data/08_Peak_Analysis/pretty_gene_names.csv` — gene symbols used in the figures (shipped).
- `input_data/08_Peak_Analysis/group_assignments_WOB_EXCLUDED.txt` — lineage of each fish (TP, BP1–BP4) from the
  ancestral character reconstruction (Fig. 2), used for the per-lineage pseudo-R² (shipped).

## Steps
Run from the repository root after `source config.env`.
1. `bash code/08_Peak_Analysis/01_get_SNP_HITS.sh PKINGS_ALL_WOB_EXCLUDED` — genotypes (0/1/2) of the TOP250 peak SNPs
   for all 305 samples (`.bin.pos.genos`). The published matrix used the 645 SNPs of the extended-window run of
   `gwas_report.Rmd`.
2. `bash code/08_Peak_Analysis/02_Render_Peak_Report.sh PKINGS_ALL_WOB_EXCLUDED 0.15` — render
   `peak_analysis_new.Rmd`: thresholds, LD-coloured regional association plots, overlapping genes and genotype
   matrices (Fig. 3B).
3. `bash code/08_Peak_Analysis/03_Render_PVE_report.sh PKINGS_ALL_WOB_EXCLUDED 0.15` — render `PVE_new.Rmd`:
   McFadden's pseudo-R² of each index SNP and SNPs in LD (r² > 0.6) with it, per BP lineage against all TP fish
   (Fig. 3C).

## Outputs
- `output_data/08_Peak_Analysis/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED.bin.pos.genos` — genotype matrix.
- `output_data/08_Peak_Analysis/PKINGS_ALL_WOB_EXCLUDED/*.peak_analysis_*.svg`, `*.peak_report.docx`,
  `*_Peak_Analysis.report.nb.html` — Fig. 3B panels and peak report.
- `output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_GENES_IN_PEAKS.txt` — genes in or near peaks.
- `output_data/08_Peak_Analysis/PKINGS_ALL_WOB_EXCLUDED/mcfaddenrsquared.svg`, `mcfaddenrsquared.dataset.Rdata`,
  `*_PVE_Report.report.nb.html` — Fig. 3C.
