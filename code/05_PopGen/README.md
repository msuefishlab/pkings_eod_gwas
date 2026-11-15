# 05_PopGen: windowed nucleotide diversity, differentiation and Tajima's D

Produces: Fig. S1 (π, Tajima's D and pairwise FST per population); the π windows are also read by the
reference-bias check (03_VCF_QC_Filter) and by later modules.

## Inputs
- `output_data/03_QC/${allsites_vcf}` — filtered biallelic SNPs plus invariant sites; `output_data/03_QC/${snps_only_vcf}`
  — filtered biallelic SNPs (03_VCF_QC_Filter).
- `input_data/05_PopGen/popsfile.txt` — sample ID and population (species, or creek for *P. kingsleyae*), tab-separated.
- `input_data/00_Reference_Genome/${reference}.fai` — one Tajima's D array task per sequence.
- `input_data/01_Terra/data_model/fish_data_2024a.txt` — sample sheet (`Rscript code/00_common/make_sample_sheets.R`).

## Steps
Run from the repository root after `source config.env`; edit the `#SBATCH` headers for your cluster.
1. `bash code/05_PopGen/01_run_pixy.sh` — pixy (`pixy.sb`): π, dxy and Weir and Cockerham's FST in non-overlapping
   5-kb windows.
2. `bash code/05_PopGen/02_run_tajd.sh` — vcftools Tajima's D in 5-kb windows for each population, one sequence per
   array task (`tajd.sb`, `split_popfile.sh`).
3. `bash code/05_PopGen/03_merge_tajd.sh`, then `bash code/05_PopGen/04_merge_pop_tajd.sh` — merge the chunks into one
   table per population and then one table for all populations.
4. `bash code/05_PopGen/05_run_popgen_stats_report.sh` — renders `overall_popgen_stats.Rmd`
   (`render_popgen_report.sb`, `render_popgen_report.R`), which writes Fig. S1 as `popstats.svg`. Its last chunks
   compare FST with river distance between collection sites; they need the R packages riverdist and sf and GIS
   layers (`input_data/XX_Mapping/`) that are not distributed, and they are not part of Fig. S1.

## Outputs
- `output_data/05_PopGen/${popgen_prefix}_pi.txt`, `_dxy.txt`, `_fst.txt` — pixy windows
- `output_data/05_PopGen/<POP>.tajd.stats.txt`, `${popgen_prefix}_tajd.txt` — Tajima's D windows
- `output_data/05_PopGen/popstats.svg`, `${popgen_prefix}.report.report.html` — Fig. S1 and report
