# 10_Introgression: ABBA-BABA tests for shared BP alleles at the GWAS peaks

Produces: Fig. S5 and Table S4 (`DsuiteAnalysis.Rmd`), Fig. S6 (`DsuiteInvestigateAnalysis.Rmd`), the Peak 6 CURV3
test (`06_peak6_curv3/`), the BP2–BP3 relative node depth (`05_bp2_bp3_dxy_at_peaks.R`) and the statistics in
Supplementary Data S7 (ESM, Analysis of Introgression).

## Inputs
- `input_data/10_Introgression/SETS.txt` — sample → group (BP1, TP1, BP2, TP2, BP3, TP_OUT, BP_ANC, Outgroup;
  209 fish) (shipped).
- `input_data/10_Introgression/pops.tree` — `((((((BP1,TP1),((BP2,TP2),BP3)),TP_OUT),(BP_ANC)),OUTGROUP));` (shipped).
- `input_data/10_Introgression/peaks.bed` — the six extended peak windows (−log10 p > 3 merged within 75 kb; the
  extended-window `PEAKS.bed` of `06_Association/gwas_report.Rmd`) (shipped; also read by 11, 12 and 15).
- `input_data/10_Introgression/test_trios.txt` — Dinvestigate trios, one tab-separated `P1 P2 P3` per line (not
  shipped). The published windows cover nine trios: BP1 TP1 BP2, BP1 TP1 BP3, BP2 TP2 BP3, BP1 TP_OUT BP_ANC,
  BP2 TP_OUT BP_ANC, BP3 TP_OUT BP_ANC, BP2 TP1 BP1, BP3 TP2 BP1 and BP3 TP2 BP2. Two of them (BP2 TP1 BP1 and
  BP3 TP2 BP2) have P1 and P2 that are not sisters and are not used by `DsuiteInvestigateAnalysis.Rmd`.
- `output_data/03_QC/${snps_only_vcf}` — unphased biallelic SNPs (03_VCF_QC_Filter).
- `output_data/05_PopGen/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites_dxy.txt` — pixy 5-kb dxy (05_PopGen), for 05.
- Container: `${root}/images/dsuite_tools.sif` (Dsuite 0.5 r58), built from `images/DsuiteDocker` (`images/README.md`).

## Steps
Run from the repository root after `source config.env`. Edit the `#SBATCH` lines of the `.sb` scripts for your cluster.
1. `bash code/10_Introgression/01_submit_dsuite.sh` — Dsuite Dtrios per chromosome, chr1–25 (`run_dsuite.sb`;
   writes `output_data/10_Introgression/chrN_*`).
2. Combine the chromosomes, as in `06_peak6_curv3/03_run_dtrios.sh`:
   `singularity exec images/dsuite_tools.sif /opt/Dsuite/Build/Dsuite DtriosCombine -o output_data/10_Introgression/whole_genome -t input_data/10_Introgression/pops.tree $(for c in {1..25}; do echo output_data/10_Introgression/chr$c; done)`
   — writes `whole_genome_combined_tree.txt`, the genome-wide D of every trio.
3. `bash code/10_Introgression/03_submit_dsuite_peaks.sh` — Dtrios within each extended peak window
   (`run_dsuite_peaks.sb`; `output_data/10_Introgression/peaks_extended/peak_N_tree.txt`).
4. `bash code/10_Introgression/04_submit_dinvestigate.sh` — Dinvestigate `-w 50,25` per chromosome for the trios in
   `test_trios.txt` (`run_dsuite_investigate.sb`). It reads bgzipped per-chromosome SNP VCFs
   `${scratch_dir}/dsuite/vcfs/chrN.vcf.gz` cut from `${snps_only_vcf}` (e.g. `bcftools view -r chrN -Oz`). Fig. S6
   uses all 25 chromosomes: use the commented `for chr in {1..25}` loop; the active loop covers only the six peak
   chromosomes.
5. Render `DsuiteAnalysis.Rmd` — one trio per hypothesis, D signed so + = sharing between the named lineages, BH q
   over the 48 peak tests, comparison with the genome-wide D (Fig. S5, Table S4).
6. Render `DsuiteInvestigateAnalysis.Rmd` — f_dM medians within each extended peak against single background windows
   and equal-size blocks of background windows (Fig. S6).
7. `Rscript code/10_Introgression/05_bp2_bp3_dxy_at_peaks.R` — RND between BP2 and BP3 at Peaks 2, 3 and 6 against
   equal-size genome blocks.
8. `code/10_Introgression/06_peak6_curv3/` — Peak 6 test with CURV3 as the third group (its own README).

Notebooks render with, e.g., `Rscript -e 'rmarkdown::render("code/10_Introgression/DsuiteAnalysis.Rmd")'`.

## Outputs
`output_data/10_Introgression/`:
- `chrN_tree.txt`, `chrN_BBAA.txt`, `chrN_Dmin.txt`, `whole_genome_combined_*.txt` — genome-wide Dtrios.
- `peaks_extended/peak_N_tree.txt` — Dtrios per peak.
- `*_localFstats_chrN_50_25.txt` — Dinvestigate windows.
- `peak_dtrios_summary.tsv` (Table S4), `dsuite_six_hypotheses_by_peak.csv`, `figures/FigS5_peak_D.{pdf,png}`.
- `peak_fdm_summary.tsv`, `figures/FigS6_fdM.{pdf,png}`.
- `bp2_bp3_dxy_peaks.tsv` — RND and dxy at the shared BP2/BP3 peaks.
- `peak6_curv3/` — the CURV3 test.

The published per-peak values are Supplementary Data S7 (`SupplementalDataS7_peak_ABBA_BABA.tsv`).

## Glossary
- **D (ABBA–BABA)** — excess shared derived alleles between P3 and one of the sisters P1, P2; D ≈ 0 under incomplete
  lineage sorting alone.
- **Z, p** — block-jackknife significance of D.
- **f4-ratio** — admixture-proportion estimate for the trio.
- **f_d, f_dM, d_f** — window statistics from Dinvestigate; f_dM is symmetric around zero.
- **BBAA / ABBA / BABA** — counts of site patterns; low counts make D unstable.
