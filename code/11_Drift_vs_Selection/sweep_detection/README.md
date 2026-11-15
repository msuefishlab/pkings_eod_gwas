# sweep_detection: calibrated tests for selective sweeps at the six GWAS peaks

Produces: Supplementary Fig. S7 and Table S5 (with `../folded_clr/` for the folded-CLR row and the six-test
correction); ESM "Tests for Selection at GWAS Peaks". `evidence/DESIGN.md` describes the statistics, the search
operator, the matched null, the decision rule and the sensitivity rows.

Three statistics per peak, each in the lineage(s) where its BP allele arose against the matched TP reference:
Rsb (`rehh`), the SweepFinder2 CLR, and πBP/πTP (pixy). Each is reduced to one value per ±100-kb search region
and ranked against ~4,400 non-peak regions reduced the same way and matched on SNP density and callable fraction.

## Inputs
- `output_data/07_Phasing/${phased_vcf}` — SHAPEIT4-phased SNP VCF, from `code/07_Phasing`
- `output_data/03_QC/${allsites_vcf}` (+ `.tbi`) — filtered SNPs + invariant sites, from `code/03_VCF_QC_Filter`
  (read by `pixy_lineages/`)
- `output_data/05_PopGen/${popgen_prefix}_pi.txt` — per-population pixy run from `code/05_PopGen` (read only by
  `pixy_lineages/validate_pixy.R`)
- `input_data/11_Drift_vs_Selection/sweep_detection/` — sample lists: `BP1_TP1.txt`, `BP2_TP2.txt`, `BP3_TP2.txt`
  (rehh pairs), `split_up/{BP1,TP1,BP2,TP2,BP3}.txt` (lineage groups), `popsfile_lineages.txt` (pixy populations)
- `input_data/10_Introgression/peaks.bed` — extended GWAS peaks and index SNPs
- `input_data/00_Reference_Genome/${reference}.fai` — reference index (external, see the root README)
- `images/{gwas_tools,rehh,sweepfinder2,dadi_tools}.sif` — built from `images/` (see `images/README.md`)

`bash code/11_Drift_vs_Selection/sweep_detection/00_fetch_inputs.sb` checks that these exist and links the all-sites
VCF to `${sweep_scratch}/03_QC/` for `../simulation/`.

## Steps
Run from the repository root after `source config.env`. Submit `.sb` files with `sbatch --export=root=$(pwd)`;
edit the `#SBATCH` headers for your cluster. Container steps need a node with singularity.

1. `bash code/11_Drift_vs_Selection/sweep_detection/00_fetch_inputs.sb` — input check.
2. `sbatch code/11_Drift_vs_Selection/sweep_detection/ehh/00_prepare_vcfs.sb` — polarize the phased VCF from the
   PSZA_ outgroup (`ehh/01_polarize_vcf.sh`), split it into the pair and lineage sets
   (`ehh/02_split_by_sample_sets.sh`) and by chromosome (`ehh/03_split_by_chromosome.sh`). Where a job is limited
   to 4 h, run `ehh/01_polarize_vcf.sh` inside `gwas_tools.sif`, then the arrays `ehh/01a_split_sets_array.sb`
   and `ehh/01b_split_by_chrom_array.sb` (same bcftools calls, one task per set).
3. These are independent of one another once step 2 has finished:
   - `sbatch .../ehh/14_ehh_scan_array.sb` — rehh inES and Rsb, 3 pairs × 25 chromosomes, via
     `ehh/scan_ehh_integrated.R`; saves the full `scan_hh` output as RDS under `${sweep_scratch}`. Run it a second
     time with `--export=root=$(pwd),SUBSAMPLE=16` for the down-sampled row.
   - `sbatch .../sweepfinder2/10_sf2_convert.sb`, then `11_sf2_background.sb`, then `12_sf2_scan.sb` — CLR on a
     10-kb grid for the five lineages (`vcf_to_sweepfinder2.py` converts the VCFs).
   - `sbatch .../pixy_lineages/run_pixy_lineages.sb` — π, dxy and Fst in 5-kb windows for the five lineage groups;
     then `Rscript .../pixy_lineages/validate_pixy.R` checks BP3 against population BIR of the 05_PopGen run.
   - `sbatch .../ehh/17_xpnsl_array.sb`, then `python3 .../ehh/17b_xpnsl_merge.py <raw_dir> <out_dir>` with
     `output_data/11_Drift_vs_Selection/sweep_detection/ehh/xpnsl_raw` and `.../ehh/xpnsl` — XP-nSL for its
     sensitivity row.
4. `bash code/11_Drift_vs_Selection/sweep_detection/evidence/00_run_all_downstream.sh [--wait]` — merges the CLR
   (`sweepfinder2/13_sf2_merge.sh`) and rehh scans (`ehh/05_merge_ehh_results.sh`), computes iHS, XP-EHH and a
   genome-wide Rsb (`ehh/10_genome_stats.R`), then runs `evidence/01_index_snp_facts.R` (BP-allele state and the
   applicability map), `02b_region_callable.R` (matching covariates), `02_region_stats.R` (the search operator on
   peaks and pseudo-peaks), `03_evidence_matrix.R` (percentiles, calls, empirical p, calibration) and
   `04_run_sensitivity.sh` (± 50-kb region, 40-kb window, unmatched and differentiation-matched nulls; the
   genome-wide-Rsb and n = 16 rows once their inputs exist).
5. Further sensitivity rows of Table S5:
   - `sbatch .../evidence/06_recomb_null.sb` — recombination-proxy-matched null (`02c_region_ines.R` + 03).
   - `sbatch .../evidence/05_n16_downstream.sb` — the n = 16 row and its `_rsbgw` control, after the SUBSAMPLE=16 scan.
   - `sbatch --array=0-7 .../evidence/08_ext_and_xpnsl.sb` — extended-peak-interval row (tasks 0–5), XP-nSL in
     place of Rsb (task 6) and a regression check of the headline matrix (task 7).
   - Folded-spectrum CLR: `../folded_clr/` (`MODE=folded` on steps 3's SweepFinder2 scripts writes the same inputs).
   - Centred 40-kb window (Fig. S7 legend), run inside `gwas_tools.sif`:
     `Rscript .../evidence/02_region_stats.R --half_width 20000 --window 40000 --suffix _centred`,
     `Rscript .../evidence/02b_region_callable.R --half_width 20000 --suffix _centred`,
     `Rscript .../evidence/03_evidence_matrix.R --suffix _centred`.
6. `Rscript code/11_Drift_vs_Selection/sweep_detection/evidence/04_evidence_figure.R` — Fig. S7 (runs locally;
   needs data.table, ggplot2, patchwork, ragg).

Notes:
- `ehh/15_verify_against_published.R` compares `ehh/merged/` with a reference scan in `ehh/merged_published/`.
  `00_run_all_downstream.sh` skips it when no reference is present; `05_n16_downstream.sb` requires it, so copy
  `ehh/merged/` to `ehh/merged_published/` after the full-sample merge and before the n = 16 run. Do not write new
  scans into `merged_published/`.
- `ehh/12_nsl.py` computes nSL, a secondary statistic `02_region_stats.R` reads if present; it is not in the
  decision rule or Table S5.
- inES standardized to mean 0 and SD 1 is written as `inES_z`; real iHS comes from `10_genome_stats.R`.
- pixy must stay at 1.2.10.beta2 (`images/Dockerfile`).
- Large intermediates (scan RDS files) live under `${sweep_scratch}`.

## Outputs
In `output_data/11_Drift_vs_Selection/sweep_detection/`:
- `MIC4273_...polarized.vcf.gz`, `AA.tsv.gz`, `outgroup_counts.tsv`, `lists/`, `vcfs/`, `vcfs/by_group/`, `by_chrom/` — polarized VCFs
- `ehh/{scans,merged,ihs,xpehh,rsb_check,xpnsl}/` — rehh and XP-nSL statistics
- `sweepfinder2/{sfs_input_unfolded,background,clr_unfolded,merged}/` — CLR
- `pixy_lineages/` — lineage-group π, dxy, Fst
- `evidence/` — `index_snp_facts.tsv`, `applicability_map.tsv`, `region_stats*.tsv.gz`, `region_callable*.tsv`,
  `evidence_matrix*.tsv` (Table S5; one file per sensitivity row), `calibration*.tsv`,
  `statistic_correlations*.tsv`, `FigS7_evidence_matrix.{pdf,png}`
