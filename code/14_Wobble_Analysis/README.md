# 14_Wobble_Analysis: mixed-waveform ("wobble") fish at the GWAS peaks

Produces: Fig. S4 (phased haplotypes at the six peaks), Table S2 (peak genotypes of mixed fish), and the Supplementary
Note analyses of mixed fish: peak haplotype states, sharing of BP haplotypes outside their home populations, ancestry
of mixed fish, the Bambomo P0 scans with manual and classifier labels, the recurrent wobble leads, the power
simulation, and the Apassa ancestry distribution (`apassa_ancestry/`, exploratory; its own README). The code calls the
paper's "Mixed" fish "Wobble".

## Inputs
- `input_data/01_Terra/data_model/fish_data_2023.txt`, `fish_data_2024a.txt` — sample sheets
  (`Rscript code/00_common/make_sample_sheets.R`)
- `output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS.txt`, `_PEAKS.bed` — peak
  SNPs and peaks of the main scan (`06_Association/06_Render_GWAS_report.sh PKINGS_ALL_WOB_EXCLUDED`)
- `output_data/07_Phasing/${phased_vcf}` — phased SNPs (`07_Phasing`)
- `output_data/08_Peak_Analysis/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS_TOP250.vcf.gz` — unphased
  genotypes (GT:AD:DP:GQ) at the Fig. S4 SNP set, the 645 extended (−log10 p > 3) peak SNPs listed in
  `input_data/15_TFBS_Scan/PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS.txt`; make it as in `code/15_TFBS_Scan/README.md`, step 1
- `output_data/03_QC/...FILTERED_BIALLELIC_SNPS_ONLY.renamed.vcf.gz` or `...FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz`
  — unphased call set for the optional raw-call check in `16` (`03_VCF_QC_Filter`)
- `output_data/06_Association/BAM_ONLY_BP0_TPWOB1/` — genotypes, SNP map, kinship, `.order` and
  `.prune.for_gemma.geno` of the 67 Bambomo fish (`06_Association` steps 1–2 for `BAM_ONLY_BP0_TPWOB1`)
- `output_data/02_Phenotyping/p0_statistical_table.csv`, `measurement_data.csv` — P0 classifier calls and EOD
  measurements (`02_Phenotyping`; they need the raw EOD recordings, which are not distributed)
- `input_data/14_Wobble_Analysis/BAM_P0/` — P0 phenotype columns for the Bambomo scans (`BAM_P0.pheno`,
  `BAM_P0.pheno_columns.tsv`, `BAM_P0.samples.tsv`), written by `07_make_p0_phenos.R` from the classifier table and
  shipped here so that the scans and the power simulation (steps 10, 11, 13, 14) run without the raw recordings
  (steps 3 and 12 also read the `02_Phenotyping` tables above)
- `output_data/13_Structure/kings_only/` — ADMIXTURE K = 8 and K = 10 (`13_Structure`)
- `input_data/04_Phylogenomics/pkings_iqtree_genes_astral_gcf_scf_tree.cf.tree.ultra.time.nex` — population tree (time-calibrated; shipped because it was written interactively from the 04 Paramormyrops notebook)
  (`04_Phylogenomics`)
- `input_data/00_Reference_Genome/paramormyrops_ncbi_0.4.renamed.sorted.gff` — the renamed NCBI annotation sorted by
  chromosome and start (read by `10` for the genes in each LD block)

## Steps
Run from the repository root after `source config.env`. R scripts need tidyverse and data.table (18 also ape, ggtree,
patchwork); scripts 00, 10, 16, 17 and 19 call `bcftools` from the PATH. `08` and `power_sim.sb` submit SLURM jobs; edit
the `#SBATCH` headers of `p0_gemma.sb` and `power_sim.sb` for your cluster.

Peak haplotypes, Fig. S4 and Table S2:
1. `bash code/14_Wobble_Analysis/00_extract_phased_peak_genotypes.sh` — phased genotypes at the peak SNPs, focal fish.
2. `Rscript code/14_Wobble_Analysis/01_call_peak_haplotypes.R` — BP/TP haplotype call per chromosome copy and peak.
3. `Rscript -e 'rmarkdown::render("code/14_Wobble_Analysis/wobble_multilocus_haplotype_analysis.Rmd", output_dir = "output_data/14_Wobble_Analysis")'`
   — multilocus combinations of the six peak states in Bambomo BP vs mixed fish.
4. `Rscript code/14_Wobble_Analysis/15_mixed_peak_genotype_table.R` — Table S2, heterozygosity and ADMIXTURE ancestry
   of mixed fish.
5. `Rscript code/14_Wobble_Analysis/16_nonfocal_peak_haplotypes.R` — the same haplotype scoring for every fish in the
   phased VCF.
6. `Rscript code/14_Wobble_Analysis/17_peak_hap_match_audit.R` — grades BP-class haplotypes outside their home
   populations; scores them on the Fig. S4 SNP set.
7. `Rscript code/14_Wobble_Analysis/19_haplotype_identity_test.R` — same haplotype vs a similar one.
8. `Rscript code/14_Wobble_Analysis/18_figS4_replacement.R` — Fig. S4.

Bambomo P0 scans, recurrent leads and power:
9. `Rscript code/14_Wobble_Analysis/07_make_p0_phenos.R` — rewrites `input_data/14_Wobble_Analysis/BAM_P0/` (needs
   the classifier table; skip it to use the shipped files).
10. `bash code/14_Wobble_Analysis/08_run_p0_gemma.sh 0`, then `bash code/14_Wobble_Analysis/08_run_p0_gemma.sh "1-100%25"`
    — GEMMA Wald scan of each phenotype column, then 100 permutations each (`p0_gemma.sb`).
11. `bash code/14_Wobble_Analysis/09_reduce_p0_gemma.sh`, then `Rscript code/14_Wobble_Analysis/09_collect_p0_gemma.R`
    — permutation thresholds, peaks and landmark p values per scan.
12. `Rscript -e 'rmarkdown::render("code/14_Wobble_Analysis/10_recurrent_leads.Rmd", output_dir = "output_data/14_Wobble_Analysis/RECURRENT_LEADS")'`
    — re-tests the five recurrent wobble leads in Bambomo and Apassa.
13. `Rscript code/14_Wobble_Analysis/11_power_sim_make.R output_data/06_Association/BAM_ONLY_BP0_TPWOB1/BAM_ONLY_BP0_TPWOB1.prune.for_gemma.geno`
    — planted-SNP simulations (prints the number of sims, N).
14. `sbatch --array=1-<ceil(N/CHUNK)> --export=root=${root},CHUNK=<CHUNK> code/14_Wobble_Analysis/power_sim.sb`, then
    `Rscript code/14_Wobble_Analysis/11_power_sim_summarize.R` — GEMMA test of each planted SNP; power against the
    permutation thresholds from step 11.

Apassa ancestry: see `apassa_ancestry/README.md`.

## Outputs
- `output_data/14_Wobble_Analysis/PKINGS_WOB_MULTILOCUS/` — `peak_haplotype_states.csv`, `peak_haplotype_allele_key.csv`,
  `multilocus_peak_haplotype_wide.csv`, `mixed_peak_state_counts.tsv` (Table S2), `mixed_admixture.tsv`;
  `nonfocal_haplotypes/`; `peak_hap_zoom/` (`identity_test_figS4snps.tsv`, `figS4_replacement.png/.pdf` = Fig. S4)
- `output_data/14_Wobble_Analysis/BAM_P0/` — per-scan GEMMA results, `SUMMARY/`, `POWER_SIM/`
- `output_data/14_Wobble_Analysis/RECURRENT_LEADS/` — lead tables and report
- `output_data/14_Wobble_Analysis/apassa_ancestry/`
