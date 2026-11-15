# local_ancestry: Bambomo ancestry along the chromosomes of the 2019 Apassa mixed fish

Produces: the exploratory local-ancestry follow-up of the Apassa Bambomo-cluster component (Supplementary Note,
"Ancestry"): per-fish Bambomo-like fraction and tracts, the simulation calibration, agreement with the
genotype projection score, and the test of enrichment at the six GWAS peaks. Scripts 09–14 are post hoc
descriptive follow-ups.

## Inputs
- `input_data/01_Terra/data_model/fish_data_2023.txt` — sample sheet; generate it first with
  `Rscript code/00_common/make_sample_sheets.R` (from Supplementary Data S1).
- `output_data/11_Drift_vs_Selection/sweep_detection/{APA,BAM,BEN}.polarized.vcf.gz` (+ `.csi`) — per-population
  subsets of the polarized SHAPEIT4 VCF written by `sweep_detection/ehh/01_polarize_vcf.sh` (via
  `sweep_detection/ehh/00_prepare_vcfs.sb`). Populations are taken by sample-name prefix (`APA_`, `BAM_`, `BEN_`;
  so `APA_223`, a Bambomo fish, sits in the APA file, and 01 regroups it from the metadata). No script here writes
  them; make them with the same bcftools steps as `02_split_by_sample_sets.sh`:
  ```bash
  pol=output_data/11_Drift_vs_Selection/sweep_detection/${phased_vcf%.vcf.gz}.polarized.vcf.gz
  lists=output_data/11_Drift_vs_Selection/local_ancestry/pop_lists; mkdir -p ${lists}
  for p in APA BAM BEN; do bcftools query -l ${pol} | grep "^${p}_" > ${lists}/${p}.txt; done
  bash code/11_Drift_vs_Selection/sweep_detection/ehh/02_split_by_sample_sets.sh ${lists} .
  ```
- `output_data/11_Drift_vs_Selection/simulation/26_ancestry_check.tsv` — projection score q per fish, from
  `simulation/26c_ancestry_check.py` (01 uses it to define the clean Apassa TP panel).
- `output_data/11_Drift_vs_Selection/simulation/26_temporal/` (`samples.txt`, `geno/`) — unphased genotype tables
  from `simulation/26a_extract_temporal.sh` (read by 11 only).
- `output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz` —
  unphased filtered callset from `03_VCF_QC_Filter` (`05_combine_filtered.sh`, then `08_rename_vcfs.sh`); passed to 06.
- `output_data/07_Phasing/${phased_vcf}` — phased VCF from `07_Phasing` (read by 14).
- `output_data/08_Peak_Analysis/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS_TOP250.vcf.gz` — from
  `08_Peak_Analysis/01_get_SNP_HITS.sh` (read by 14).
- `input_data/10_Introgression/peaks.bed` (extended peak intervals) and
  `input_data/00_Reference_Genome/${reference}.fai`.

## Steps
Run from the repository root after `source config.env`, with python3 (numpy, scipy) and bcftools on the path
(locally or inside `gwas_tools.sif`); about 15 minutes in total. `OUT=output_data/11_Drift_vs_Selection/local_ancestry`.
1. `python3 code/11_Drift_vs_Selection/local_ancestry/01_panels.py $OUT` — groups from the metadata, reference
   panels, clean Apassa TP (q < 0.05), hold-outs (seed 20261103).
2. `bash code/11_Drift_vs_Selection/local_ancestry/02_extract_candidates.sh $OUT 10` — merges the APA/BAM/BEN phased
   VCFs and keeps sites with |Δp| ≥ 0.5 under any panel pair.
3. `python3 .../03_markers.py $OUT` — ancestry-informative markers chosen against an independent TP sample, 10-kb
   thinning, r² pruning, empirical-Bayes emission priors.
4. `python3 .../04_fit_hap.py $OUT <panel>` for `primary`, `full` and `ben` — phased two-state HMM per fish: grid ML
   (T, m), chromosome bootstrap, posteriors, tracts, q_m, F_long.
5. `python3 .../05_simulate.py $OUT` — simulated admixed fish from held-out haplotypes on a T × m grid, fitted the
   same way; leave-one-out matching and power.
6. `python3 .../06_fit_dip.py $OUT <unphased VCF>` — three-state diploid HMM on phased-VCF genotypes and on the
   unphased calls (DP ≥ 8, GQ ≥ 20).
7. `python3 .../07_match.py $OUT` — matches real fish to simulated cells; coarse classes, site check, panel
   sensitivity, tract sharing, correlation with the 26c score.
8. `python3 .../08_peaks.py $OUT` — Bambomo-ancestry posterior at each extended peak against size-matched windows.
9. `python3 .../09_posthoc.py $OUT` — long-tract share and longest tract against simulation; per-copy peak posteriors.
10. `python3 .../10_shared_segments.py $OUT` — Bambomo-like segments shared by the 2019 wobbles, rotation null.
11. `python3 .../11_apassa_pca.py $OUT` — PCA of the 30 Apassa fish on the unphased genotype tables.
12. `python3 .../12_peak_runs.py $OUT` — Bambomo-like runs through each peak per Apassa haplotype.
13. `python3 .../13_peak5_haplotypes.py $OUT` — haplotype distances of the Bambomo-like Apassa Peak 5 copies.
14. `python3 .../14_peak5_diag_snps.py $OUT` — Peak 5 on the diagnostic SNPs of
    `code/14_Wobble_Analysis/17_peak_hap_match_audit.R`, left and right of the index SNP.

`la_lib.py` holds the HMMs, the grids, tract calling and matching (imported by 04–08).

## Outputs
- `$OUT/samples.tsv`, `groups_fill.tsv`, `cand/`, `markers_*.npz`, `markers_summary.tsv` — panels and markers
- `$OUT/hap_<panel>/`, `dip_*/`, `unphased_calls/`, `sim/` — fits, posteriors, tracts and simulations
- `$OUT/results/` — `per_fish.tsv`, `summary.txt`, `peaks.tsv` and the post hoc tables of steps 9–14
