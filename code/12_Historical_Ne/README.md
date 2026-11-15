# 12_Historical_Ne: historical effective population size with SMC++

Produces: Fig. S9 and the ESM Methods section "Historical Effective Population Size" (SMC++ 1.15.4 Ne(t) for four BP
and four TP populations, and the simulation check of what SMC++ can recover).

## Populations and samples

| Lineage | Population (label) | Fish used by `vcf2smc` |
|---|---|---|
| BP | Bambomo (`BAM`) | 20 of 66 |
| BP | Bikagala (`BIK`) | 20 of 30 |
| BP | Bavavela (`BAVA`) | 20 |
| BP | Biroundou (`BIR`) | 15 |
| TP | Apassa, 2009 sample (`APA_2009`, labelled `APA` inside SMC++) | 20 |
| TP | Bengue (`BEN`) | 17 |
| TP | Douvalou (`DOV`) | 16 |
| TP | Douengi (`DOG`) | 15 |

Population membership is from `input_data/05_PopGen/popsfile.txt`. Excluded: the ten Apassa fish sampled in 2019, one
fish of each second-degree relative pair (`BAM_175`, `BIR_336`) and the Bavavela relative cluster (`BAVA_6626`,
`BAVA_6627`, `BAVA_6890`). Each list holds at most 20 fish: Bambomo and Bikagala keep five distinguished individuals
plus 15 others drawn at random (seed 42). The first five fish of every list are the distinguished individuals (DIs);
for Bambomo these are five mixed-phenotype fish sampled in 2019.

## Inputs
- All-sites VCF `output_data/03_QC/${allsites_vcf}` (biallelic SNPs + invariant sites; from `code/03_VCF_QC_Filter`),
  split into `${scratch_dir}/smcpp/data/vcf_by_chr/chr{1..25}.vcf.gz` (+ `.tbi`) with `bcftools view -r chrN`.
- `input_data/00_Reference_Genome/${reference}.fai` (chromosome lengths).
- `input_data/10_Introgression/peaks.bed` (the six GWAS peaks; column 1 has no `chr` prefix).
- SMC++ container `images/smcpp.sif`: `singularity pull images/smcpp.sif docker://terhorst/smcpp:latest`
  (the run used SMC++ 1.15.4; check with `singularity exec images/smcpp.sif smc++ --version`).
- `input_data/12_Historical_Ne/samples/<POP>_sub20.txt` — the `vcf2smc` lists above (DIs first), read by
  `orchestrate.sh`, `w_vcf2smc.sh` and `w_joint.sh`; `samples/<POP>.txt` — the full lists after exclusions, read by
  `generate_awk.py` / `build_masks.py`.
- `input_data/12_Historical_Ne/vcf2smc_tasks_sub20_lpt.txt` — 1,000 lines `POP CHR DI DI_INDEX` (8 populations x 25
  chromosomes x 5 DIs, largest chromosomes first), read by `orchestrate.sh`.
- `${scratch_dir}/smcpp/masks/<POP>.mask.bed.gz` (+ `.tbi`; working tree overridable with `SMCPP_WORK`), the original
  per-population masks, built before step 1 and not produced by a script here: the union of sites with population
  call rate < 0.8 (`compute_callrate.awk` run on `bcftools query -f '%CHROM\t%POS[\t%GT]\n'` per chromosome with
  `-v outdir=<dir> -v chrom_name=<chr>`; `build_masks.py` is the same rule in one pass), sites with FILTER other than
  PASS, the chr4 inversion-like region chr4:40,880,000-42,770,000, and each GWAS peak +/- 1 Mb
  (`expand_peaks.py peaks.bed <fai>`). The shipped `compute_callrate.awk` is the one that ran; it was generated
  before the relative exclusions in Bavavela and Biroundou, so it counts 22 and 16 fish there. Regenerating it with
  `generate_awk.py` from `samples/<POP>.txt` gives the 20- and 15-fish lists.
- For the simulation check: a Python with msprime 1.4.4, numpy and matplotlib (`MSPPY`), and HTSlib (bgzip, tabix).

## Steps
Run from the repository root after `source config.env`, on a SLURM cluster (edit the `#SBATCH` headers and the
`--account`/`--partition` arguments inside `orchestrate.sh` for your cluster).
1. `mkdir -p output_data/slurm_logs/12_Historical_Ne && sbatch code/12_Historical_Ne/scripts/build_masks_v2.sb` —
   finds every position absent from the all-sites VCF (14.15% of chr1-25; `vcf2smc` would read them as homozygous
   reference) and adds it to the original masks: `audit/masks_v2/<POP>.mask.bed.gz`, plus `callability_by_chr.tsv`
   and `mask_comparison.tsv`.
2. `OUTD=$root/output_data/12_Historical_Ne/results sbatch code/12_Historical_Ne/scripts/run_smcpp.sb` — one
   controlling job (`orchestrate.sh`) that submits and waits on SLURM arrays (about 1,400 CPU-hours):
   `w_vcf2smc.sh` (1,000 single-population `vcf2smc` tasks on masks_v2), `w_estimate.sh` (`smc++ estimate --cores 4
   --polarization-error 0.5`, i.e. folded, at mu = 1e-8 and 4e-9), `make_bootstraps.py` + `w_estimate.sh` (10
   chromosome bootstraps per population and mu), `w_joint.sh` + `w_split.sh` (two-population split times for six
   pairs; computed but not reported in the paper), then `collate.py`. `DRYRUN=1 bash
   code/12_Historical_Ne/scripts/orchestrate.sh` checks the inputs and prints task counts without submitting. The
   job is resumable: resubmit it if it stops; finished tasks are skipped.
3. `python3 code/12_Historical_Ne/scripts/plot_Ne_figure_esm.py output_data/12_Historical_Ne/results/Ne_all_pops_both_mu.tsv output_data/12_Historical_Ne/esm_fig_s9/Fig_S9_Ne`
   — Fig. S9 (mu = 1e-8 arm).
4. `sbatch -o output_data/slurm_logs/12_Historical_Ne/validation_%j.log code/12_Historical_Ne/scripts/validation/validate_orchestrator.sh`
   — simulation check: msprime histories (constant size, constant size with gene flow, and a drop from Ne 47,250 to
   2,500 at 500, 1,000, 2,000, 5,000 or 10,000 generations; 20 diploids, 25 chromosomes, 20% of sites removed)
   passed through the same `vcf2smc`/`estimate` settings (`simulate_chrom.py`, `w_vcf2smc.sh`, `w_estimate.sh`), then
   `analyze_results.py`. This supports the Fig. S9 legend: reductions within about the last 1,000 generations are
   recovered, older ones are not timed correctly, and the minimum size comes out about two-fold too low.

## Outputs
- `output_data/12_Historical_Ne/results/Ne_all_pops_both_mu.tsv` — Ne(t) for 8 populations x 2 mu on 200
  log-spaced times (10 to 2 x 10^6 generations), with 95% bootstrap bands; `Ne_tables/<POP>_mu<arm>_Ne_table.tsv` per
  population; `rho_mu_table.tsv` (rho, theta, r/mu raw and divided by w = 100); `split_table.tsv`;
  `collation_status.txt`; the SMC++ model JSONs (`est__*`, `boot__*`, `split*__*`).
- `output_data/12_Historical_Ne/esm_fig_s9/Fig_S9_Ne.{pdf,png}` — Fig. S9.
- `${scratch_dir}/smcpp/audit/` — masks_v2 and callability tables (step 1).
- `${scratch_dir}/smcpp/results/validation/validation_summary.tsv`, `validation_figure.png` — simulation check.
