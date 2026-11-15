# Drift test: fixed differences and the neutral expectation

Produces: Table S6 (drift test at the GWAS index SNPs; `13_report_tables.py`) and the ESM Methods section "Drift
Test at GWAS Peaks": dadi two-population fits per BP–TP pair, the expected number E of SNPs as differentiated as
each index SNP, and the observed count O. It also writes the near-fixed-difference clusters used as a null by
`../../folded_clr/` and the fitted histories and block scans read by `../../simulation/` and `../model_extension/`.
What each script computes is in `METHODS.md` (the script comments cite its section numbers).

## Inputs
- `${sweep_scratch}/03_QC/${allsites_vcf}` (+ `.tbi`) — the filtered all-sites VCF from `03_VCF_QC_Filter`, staged
  by `11_Drift_vs_Selection/sweep_detection/00_fetch_inputs.sb`.
- `input_data/11_Drift_vs_Selection/demographic_modeling/` — `groups.tsv`, `group_sizes.tsv`, `pairs.tsv`,
  `chunks.tsv` (shipped; rewritten by step 0) and `pixy_reference_totals.tsv` (per-chromosome pixy π and dxy
  totals from `05_PopGen`, for the step 2 check).
- `input_data/05_PopGen/popsfile.txt`, `input_data/11_Drift_vs_Selection/sweep_detection/split_up/*.txt`,
  `input_data/10_Introgression/peaks.bed` (column 6 = index SNP).
- `input_data/00_Reference_Genome/`: the `.fai` of `${reference}` and `paramormyrops_ncbi_0.4.renamed.gff`
  (see the root README, "Reference genome").
- `output_data/11_Drift_vs_Selection/sweep_detection/evidence/index_snp_facts.tsv` — from `sweep_detection/evidence/`.
- Containers `${gwas_tools_image}` (bcftools) and `${dadi_image}` (dadi 2.4.3), built from `images/`.

## Steps
Run from the repository root on the cluster after `source config.env`, and first
`export FD_OUT=output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree`: every script writes
under `${FD_OUT}`, and the downstream modules read this directory (without it the scripts write to
`.../fixed_differences/`). Wait for each step to finish before the next unless the command chains the jobs.

| Step | Command | Produces |
|---|---|---|
| 0 | Inside `${gwas_tools_image}`: `bcftools query -l "${sweep_scratch}/03_QC/${allsites_vcf}" > vcf_samples.txt`, then `python3 code/11_Drift_vs_Selection/demographic_modeling/fixed_differences/00_prepare_inputs.py --vcf_samples vcf_samples.txt` | the input tables above; `${FD_OUT}/masks/{exons,exclude}.bed` |
| 1 | `sbatch --export=root=$(pwd),FD_OUT=${FD_OUT} --array=0-102%40 code/11_Drift_vs_Selection/demographic_modeling/fixed_differences/01_site_table.sb` | `site_table/` per 10-Mb chunk |
| 2–6 | `bash code/11_Drift_vs_Selection/demographic_modeling/fixed_differences/run_to_spectra.sh` | site-table check (02; later jobs run only if it passes), projections (04), call deserts (05), `sfs/<pair>.pkl`, `tail/index_units.tsv` (06) |
| 6b | `06b_block_scan.py` inside `${dadi_image}` (e.g. `sbatch --wrap "singularity exec --bind ${scratch_dir} ${dadi_image} python3 code/11_Drift_vs_Selection/demographic_modeling/fixed_differences/06b_block_scan.py"`) | `sfs/<pair>.block_scan.tsv` (E_hetNe weights) |
| 7 | `bash .../run_fits.sh fit`, then `run_fits.sh refine 2`, `nested 1`, `refine 3`, `nested 2` (each after the previous finishes) | `fits/<pair>.<model>.data.seed*.tsv` |
| 8 | `bash .../run_fits.sh collect` | `fits/best_fits.tsv`, `fits/lrt.tsv` |
| 9 | `bash .../run_fits.sh checks` (10_checks.py) | `checks/{spectra,model}_checks.tsv` |
| 10 | `bash .../run_fits.sh tail`, then `run_fits.sh boot` (bootstrap and Poisson refits), then `run_fits.sh tail` again | final `tail/*.tsv` |
| 11 | `11_clusters.sb` as an array over chunks, then once without `--array` (commands in its header) | cluster BEDs, `fixed_diff_cluster_counts.tsv`, `phaseB_gate.tsv`, `parallel_count.tsv` |
| 12 | `Rscript .../12_tail_figure.R ${FD_OUT}/tail <out_prefix>` | tail figure |
| 13 | `python3 .../13_report_tables.py ${FD_OUT}/tail > tables.md` | Table S6 and report tables |

`.../` is `code/11_Drift_vs_Selection/demographic_modeling/fixed_differences/`. `run_fits.sh` and `run_to_spectra.sh`
submit their own SLURM jobs; edit the `#SBATCH` headers and `sbatch` resource flags for your cluster.

Notes:
- Genotypes with no reads are set to missing in step 1 (GATK writes 0/0 with DP = 0 at invariant sites); without
  this the comparison counts run 2.6–5% above pixy's.
- Do not use `output_data/11_Drift_vs_Selection/sweep_detection/AA.tsv.gz` for polarization: it came from the
  imputed phased VCF. The ancestral allele here comes from the six outgroup species in the site table.

## Outputs
- `output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree/` — `site_table/`, `checks/`,
  `sfs/`, `fits/`, `tail/` (`unit_results.tsv` holds Δ, O, E, E_robust and the tier per unit), the cluster files
  and `fixed_diff_counts.tsv`.
