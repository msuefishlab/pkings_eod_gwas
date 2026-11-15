# 13_Structure: ADMIXTURE ancestry of the P. kingsleyae samples

Produces: Fig. S8 (ADMIXTURE ancestry proportions and cross-validation error); the K = 8 and K = 10 runs are also read
by `14_Wobble_Analysis/15_mixed_peak_genotype_table.R` (ancestry of mixed-phenotype fish, Supplementary Note).

## Inputs
- `output_data/03_QC/${snps_only_vcf}` — filtered biallelic SNPs (`03_VCF_QC_Filter`)
- `input_data/00_Reference_Genome/${chr_name_map}` — chromosome name map (see the root README, "Reference genome")
- `input_data/06_Association/PKINGS_TP1_BP2_WOB9.pheno.txt`, `.keep.txt` — the 240 P. kingsleyae (`06_Association`)
- `input_data/01_Terra/data_model/fish_data_2024a.txt` — sample sheet (`code/00_common/make_sample_sheets.R`)

## Steps
Run from the repository root after `source config.env`. `01` and `02` submit SLURM jobs; edit the `#SBATCH` header of
`prep_admixture.sb` and the `-A` account in `02_run_admixture.sh` for your cluster.
1. Make the sample list for the run `kings_only` (all 240 P. kingsleyae; PLINK wants FID and IID):
   `mkdir -p output_data/13_Structure/kings_only && awk '{print $1, $1}' input_data/06_Association/PKINGS_TP1_BP2_WOB9.keep.txt > output_data/13_Structure/kings_only/kings_only.keep.txt`
2. `bash code/13_Structure/01_run_prep_admixture.sh` — SNPs with ≤ 20% missing data, MAF ≥ 0.05 and INFO/DP ≥ 5, thinned
   to one random SNP per 20 kb (seed 42), written as a PLINK file set (`prep_admixture.sb`).
3. `bash code/13_Structure/02_run_admixture.sh kings_only 11` — subsets the PLINK file to the run's samples and runs
   ADMIXTURE with cross-validation for K = 1–11, one job per K.
4. `bash code/13_Structure/03_make_admixture_datafile.sh kings_only` — collects the Q files and CV errors.
5. `bash code/13_Structure/04_Render_ADMIXTURE_report.sh kings_only` — renders `Admixture_Results.Rmd`: CV error by K
   and the ancestry bar plot with the phenotype strip. The K drawn is set by `MAXK` in the Rmd.

## Outputs
- `output_data/13_Structure/*.for_admixture.{bed,bim,fam}` — LD-thinned SNP set, all samples
- `output_data/13_Structure/kings_only/` — per-K `.Q`, `.P`, `.out`; `cv_errors.txt`; `admixture_long_results.header.txt`;
  `kings_only_admixture.report.html`
- `output_data/13_Structure/kings_only_final_admixture.pdf` — the bar plot
