# 15_TFBS_Scan: predicted transcription-factor binding changes at GWAS peak variants

Produces: the matched-control test of directional motif changes at the three TFBS peaks (Methods, "Tests for
directional enrichment") and Supplementary Data S4–S6 (candidate variants and predicted motif changes at Peak 4,
cntnap5; Peak 5, sptbn4b; Peak 2, itih3). Data S4–S6 are `candidates/data_table_{cntnap5,sptbn4b,itih3}.tsv`.

## Inputs
- `input_data/15_TFBS_Scan/PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS.txt` (`$TFBS_SNP_METADATA`, shipped) — the 645 SNPs
  with −log10 p > 3 in the six extended peaks of the main scan (GEMMA LRT, `PKINGS_ALL_WOB_EXCLUDED`), in the column
  layout of `06_Association/gwas_report.Rmd` (`chr rs ps … allele1 allele0 af … p_lrt … log_p … peak`; `ps` is 1-based,
  allele1 is REF)
- `$TFBS_VCF` — genotypes of all 305 fish at those SNPs; made in step 1
- `input_data/15_TFBS_Scan/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt` (`$JASPAR_MOTIFS`) — JASPAR 2024
  CORE vertebrates, non-redundant, MEME format (879 motifs). Third-party file, not shipped; download it from
  <https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt>
- `input_data/15_TFBS_Scan/jaspar2024_families.tsv` — TF class and family per motif; made in step 2
- `output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED.lmm2.assoc.txt` — genome-wide GEMMA
  output, the control pool (`06_Association`)
- `input_data/10_Introgression/peaks.bed` — extended peak intervals (numeric chromosomes: `24` = `chr24`)
- `$TFBS_REF` (with `.fai`), `$TFBS_GFF` — renamed reference and NCBI annotation (`input_data/00_Reference_Genome`)
- `$phased_VCF` — phased SNPs (`07_Phasing`); used by step 8
- `output_data/03_QC/MIC4273_HAP2_NEWREF_MSU618_all_fish_all_sites.FILTERED_BIALLELIC_AND_INDEL.renamed.vcf.gz` — filtered
  SNPs and indels (`03_VCF_QC_Filter`: `05_combine_filtered.sh`, then `08_rename_vcfs.sh`); used by step 8
- `input_data/01_Terra/data_model/fish_data_2023.txt` — sample sheet (`code/00_common/make_sample_sheets.R`)

Software: MEME Suite 5.5.8 (`fimo`, `fasta-get-markov`), bcftools, samtools, python3 with numpy, pandas and scipy
(matplotlib for the optional null plots).

## Steps
Run from the repository root after `source config.env`. `$OUT` below is `output_data/15_TFBS_Scan`, where every step
writes.
1. Genotypes at the peak SNPs. This reproduces the VCF that `code/08_Peak_Analysis/01_get_SNP_HITS.sh` wrote for the
   paper from these 645 SNPs:
   ```bash
   mkdir -p output_data/15_TFBS_Scan
   tail -n +2 ${TFBS_SNP_METADATA} | cut -f1,3 > output_data/15_TFBS_Scan/peak_snps.regions.txt
   bcftools view -R output_data/15_TFBS_Scan/peak_snps.regions.txt output_data/03_QC/${snps_only_vcf} -O z -o ${TFBS_VCF}
   bcftools index ${TFBS_VCF}
   ```
   `14_Wobble_Analysis` reads the same file at `output_data/08_Peak_Analysis/PKINGS_ALL_WOB_EXCLUDED/`; copy
   `${TFBS_VCF}` and its index there.
2. `python3 code/15_TFBS_Scan/fetch_jaspar_families.py --motifs ${JASPAR_MOTIFS} --out input_data/15_TFBS_Scan/jaspar2024_families.tsv`
   — TF class and family of each motif from the JASPAR REST API (needs network access).
3. `python3 code/15_TFBS_Scan/polarize_peak_snps.py` — stage A1: BP and TP allele of each peak SNP in the TFBS windows.
4. `python3 code/15_TFBS_Scan/build_units.py` — stage A2: peak and control units and their haplotype sequences.
5. `bash code/15_TFBS_Scan/run_fimo_units.sh` — stage F: FIMO scan of all unit sequences, in 14 chunks.
6. `python3 code/15_TFBS_Scan/verify_flank_equivalence.py --units $OUT/units.tsv --hits $OUT/hits.npz --ref ${TFBS_REF} --motifs ${JASPAR_MOTIFS} --bfile $OUT/genome.bg`
   — check: the short-flank unit scans equal a full-window scan at the three peaks.
7. `python3 code/15_TFBS_Scan/matched_null.py --units $OUT/units.tsv --hits $OUT/hits.npz --families input_data/15_TFBS_Scan/jaspar2024_families.tsv --polarized $OUT/peak_snps_polarized.tsv --out_dir $OUT/results`
   — stage B: the matched-control test.
8. `python3 code/15_TFBS_Scan/candidate_variants.py` — stage C: candidate variants in LD with each index SNP and their
   predicted motif changes (Data S4–S6).

`python3 -m unittest discover -s code/15_TFBS_Scan/tests` tests `matched_null.py` on synthetic data.

## Outputs
- `$OUT/peak_snps_polarized.tsv`, `$OUT/units.tsv`, `$OUT/units.fa`, `$OUT/genome.bg`, `$OUT/hits.npz`
- `$OUT/results/` — `summary.txt`, `locus_tests.tsv`, `family_tests.tsv`, `named_group_tests.tsv`, `balance.tsv`,
  `calibration.tsv`, `fallbacks.tsv`, `per_unit_calls_{locus}.tsv`, `null_fraction_{locus|pooled}.png`
- `$OUT/candidates/` — `candidates_{locus}.tsv`, `motif_changes_{locus}.tsv`, `data_table_{locus}.tsv` (Data S4–S6),
  `candidate_haps.fa`, `summary.txt`

## Specification
The scripts implement the rules below; their docstrings give the details.

**Constants.** `W = 33`, the widest motif in `$JASPAR_MOTIFS`; `F = W − 1 = 32` bp of flank on each side, also the
largest gap between consecutive SNPs of one unit. `motif_idx` is the 0-based order of `MOTIF` lines in the motif
file. Presence call: FIMO p < 1e-4; scan threshold 1e-3 (the 1e-3 hits are used only by the strong-change rule).
Seed 20261013 throughout.

**Peak SNPs and windows.** TFBS windows, 0-based half-open (a SNP is inside if `start < ps <= end`):
sptbn4b chr17 11642886–11742887; cntnap5 chr16 10839999–11040000; itih3 chr8 15736955–15836956. Phenotypes come from
`fish_data_2023.txt` (column 1 = sample id, column 7 = `Phenotype`); only `BP` and `TP` fish are used.

**Control pool.** SNPs of the genome-wide GEMMA output with `af > 0.15` (the GWAS filter), more than 100 kb from any
extended peak in `peaks.bed`, and outside CDS (`$TFBS_GFF`, feature type `CDS` only). Exons other than CDS are allowed
in both the peaks and the pool.

**Stage A1 (`polarize_peak_snps.py` → `peak_snps_polarized.tsv`).** One row per peak SNP inside a window: `locus chr pos
snp_id ref alt log_p af_ref n_bp n_tp af_alt_bp af_alt_tp delta_af bp_allele tp_allele bp_is_ref in_cds core ref_ok`.
`af_ref` is the metadata `af` column (GEMMA `af`, which tracks the ALT allele). `af_alt_bp` and `af_alt_tp` are ALT
frequencies among called BP and TP fish; `bp_allele` is ALT if `af_alt_bp > af_alt_tp`, otherwise REF; `core` =
`log_p > 5`; `ref_ok` = REF equals the FASTA base.

**Stage A2 (`build_units.py` → `units.tsv`, `units.fa`).** A unit is one or more SNPs sorted by position; consecutive
SNPs with gap ≤ F share a unit. Peak units are clustered among the non-CDS peak SNPs of one locus, as set `extended`
(all of them) or `core` (core SNPs only). Pool units: `k = 1`, a random sample of 300,000 pool SNPs; `k ≥ 2`, up to
20,000 random chains of k consecutive pool SNPs with every gap ≤ F, for k = 2 to `KMAX = min(4, largest peak k)`. Pool
SNPs outside a unit stay reference, as non-peak SNPs do at the peak. `units.tsv` columns: `unit_idx` (0..N−1, peak units
first), `source` (`peak`/`pool`), `locus`, `set` (`extended`/`core`/`pool`), `chr`, `k`, `pos`, `ref`, `alt` (comma-
separated, ACGT), `snp_off` (1-based offsets in the unit sequence; the first is F + 1 = 33), `seq_start` (`pos_0 − F`),
`seq_len` (`pos_last − pos_0 + 2F + 1`), `gap_bin` (minimum gap binned `1-4`, `5-12`, `13-32`; `NA` for k = 1), `gc`
(GC fraction ±50 bp around the unit span, N excluded), `cpg` (per SNP, 1 if either allele makes a CpG with a reference
neighbour), `af_ref`, `snp_id`, `log_p`, and for peak rows `bp_mask` (bit j = 1 if the BP allele at SNP j is ALT).
Haplotype mask: bit j = 1 means ALT at SNP j. `units.fa` holds masks `bp_mask` and `bp_mask ^ (2^k − 1)` for each peak
unit and all 2^k masks for each pool unit, with headers `>U{unit_idx:07d}_M{mask}`; pool units whose span ±F contains
N are dropped.

**Stage F (`run_fimo_units.sh`, `parse_fimo_hits.py` → `hits.npz`).** `fimo --text --no-pgc --thresh 1e-3 --bfile
genome.bg` (0-order background from the reference). Only hits covering at least one `snp_off` of their unit are kept,
after checking that `matched_sequence` equals the unit substring. `hits.npz` holds one entry per (unit, mask, motif)
with the minimum p: `unit_idx`, `mask`, `motif_idx`, `mlog10p`, plus `motif_ids` and `motif_names`.

**Families (`jaspar2024_families.tsv`).** Columns `motif_id motif_name tf_class tf_family`; several families are
joined with `;` (the first is the primary family); an empty family becomes `Unknown`.

**Stage B (`matched_null.py` → `results/`).** Directional call per unit, motif and orientation (BP mask b, TP mask
t = b ^ (2^k − 1)): `BP-only` if the motif is present (p < 1e-4) on b and not on t, `TP-only` if the reverse. The
strong-change variant also requires no hit at p ≤ 1e-3 on the other allele. Family call: a family is present on a
haplotype if any of its motifs is. Matching of controls to peak unit u: for k = 1, stratum = (directed BP→TP class,
strand-collapsed to C>A, C>G, C>T, T>A, T>C, T>G; cpg; gc quintile; af quintile; bp_is_ref), quintile edges from the
pool, each pool SNP usable in either orientation; strata with fewer than 50 pool units drop af, then gc, then cpg
(logged). For k ≥ 2, stratum = (k capped at KMAX, gap_bin, gc quintile); per SNP j, if the pool SNP's undirected class
matches peak SNP j's, the orientation reproduces peak SNP j's directed class, otherwise peak SNP j's bp_is_ref is
copied. 10,000 draws per locus and set; draw d at each locus forms joint draw d. Two-sided p = min(1, 2 · min tail),
each tail (k + 1)/(N + 1). Statistics: (1) primary, the pooled BP-only fraction f = nBP/(nBP + nTP) over the three
loci (set `extended`); (2) per locus, f and directional events per SNP; (3) per family, nBP − nTP with BH-FDR across
families within each locus and pooled; (4) recurrence, the number of families with nBP − nTP > 0 at ≥ 2 loci;
(4b) named groups BCL6 (MA0463.3, MA0731.1) and C/EBP (MA0102.5, MA0466.4, MA0836.3, MA0837.3, MA0838.1, MA1636.2,
MA0019.2), called like families, Bonferroni across the two groups; (5) sensitivity: random sign-flip of each peak
unit's orientation (10,000 times), the strong-change rule, and set `core`; (6) calibration: 500 null draws as
pseudo-peaks must give uniform p (KS test), and the run aborts if the null SD of f is below
`0.8 · 0.5 / sqrt(mean(nBP + nTP))`; (7) power: 20 synthetic BP-only calls for one family at one locus must reach
q < 0.05; (8) balance of gc, cpg, af, class mix, k and bp_is_ref, peak vs control.

**Stage C (`candidate_variants.py` → `candidates/`).** Biallelic SNPs and indels with r² ≥ 0.8 to the index SNP and a
genotype call rate ≥ 0.8 in BP + TP fish, in the TFBS window ∪ the extended peak; phased local haplotypes ±F around
each variant; FIMO at p ≤ 1e-2. A strong change is best p < 1e-4 with a ≥ 10-fold difference between alleles. This
stage prioritises variants and computes no enrichment test.
