# 03_VCF_QC_Filter: hard filtering of the joint-genotyped call set, chromosome renaming, reference-bias check

Produces: the filtered call sets used by every later module (Methods, Variant Calling); Table S3 (reference bias).

## Inputs
- `output_data/01_Terra/${vcf_file}` and `.tbi` — raw all-sites VCF from joint genotyping (01_Terra), in the
  original scaffold names.
- `input_data/01_Terra/chrom-interval-files/`, `remaining-interval-files/` — GATK interval lists on the original
  scaffold names; regenerate them from the reference as described in `code/01_Terra/README.md`.
- `input_data/00_Reference_Genome/${original_reference}` and `${reference}` with `.fai` and `.dict`, and
  `${chr_name_map}` — see `input_data/00_Reference_Genome/README.md`.
- For the reference-bias check only: `output_data/05_PopGen/${popgen_prefix}_pi.txt` (05_PopGen),
  `output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_PEAKS.bed` (06_Association),
  `input_data/06_Association/PKINGS_ALL_WOB_EXCLUDED.pheno.txt`,
  `input_data/00_Reference_Genome/paramormyrops_ncbi_0.4.renamed.sorted.gff` and the generated sample sheet
  `input_data/01_Terra/data_model/fish_data_2023.txt`.

## Steps
Run from the repository root after `source config.env`. Scripts named `NN_*.sh` submit the matching `.sb` SLURM
job; edit the `#SBATCH` headers (account, partition, `module load` lines) for your cluster. Shards are written to
`${scratch_dir}/output_data/03_QC/`.
1. `bash code/03_VCF_QC_Filter/00_index_vcf.sh` — index the raw VCF if Terra did not return a `.tbi`.
2. `bash code/03_VCF_QC_Filter/01_submit_variant_table.sh`, then `bash code/03_VCF_QC_Filter/02_combine_interval_tables.sh`
   — per-site GATK annotations (VariantsToTable) per interval, combined into one SNP and one indel table.
3. `bash code/03_VCF_QC_Filter/03_run_qc_reports.sh` — annotation distributions (`03a_QC_Report_indels.Rmd`,
   `03b_QC_Report_snps.Rmd`) used to set the hard-filter thresholds.
4. `bash code/03_VCF_QC_Filter/04_apply_hard_filter.sh` — hard filters (SNPs: QD, SOR, FS, MQ, MQRankSum,
   ReadPosRankSum, genotype DP; indels: QD, FS, ReadPosRankSum, genotype DP and GQ), filtered genotypes set to
   no-call; biallelic SNPs, and invariant sites called in at least 50% of samples, per interval.
5. `bash code/03_VCF_QC_Filter/05_combine_filtered.sh` — concatenates the shards into
   `*.FILTERED_BIALLELIC_SNPS_ONLY.vcf.gz`, `*.FILTERED_BIALLELIC_AND_INDEL.vcf.gz` and
   `*.FILTERED_BIALLELIC_AND_INVARIANT.vcf.gz`.
6. `bash code/03_VCF_QC_Filter/06_get_gqs.sh`, then `bash code/03_VCF_QC_Filter/07_run_gq_reports.sh` — genotype
   quality report (`GQ_Report.Rmd`) on a 1% sample of sites.
7. `bash code/03_VCF_QC_Filter/08_rename_vcfs.sh` — renames scaffolds to `chr1`–`chr25` / `upsN` with
   `${chr_name_map}`, giving `${snps_only_vcf}` and `${invariant_and_snps_vcf}`.
8. Reference bias (Table S3):
   - `python3 code/03_VCF_QC_Filter/ref_bias_make_regions.py $root output_data/03_QC/ref_bias/ref_bias_regions.bed`
     — peak, high-diversity non-coding and genome-background regions;
   - `sbatch --export=root=$root code/03_VCF_QC_Filter/ref_bias_stats.sb` — per-sample depth, missing sites and
     heterozygote allele balance per region class, one chromosome per array task (`ref_bias_stats.py`);
   - `Rscript code/03_VCF_QC_Filter/ref_bias_summary.R` — BP vs TP comparisons.

`00_rename_reference.sh`, `00_rename_gff.py` and `00_rename_gtf.py` rename the reference FASTA and annotations to
the chromosome names; `input_data/00_Reference_Genome/README.md` says when they are needed.

## Outputs
- `output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_*.vcf.gz` (+ `.renamed.vcf.gz`) — call sets
- `output_data/03_QC/*.report.html` — QC and GQ reports
- `output_data/03_QC/ref_bias/ref_bias_BP_vs_TP.tsv`, `ref_bias_by_comparison.tsv`, `ref_bias_by_pop.tsv`,
  `ref_bias_peaks_by_ref_haplotype.tsv`, `ref_bias_per_sample.tsv` — Table S3
