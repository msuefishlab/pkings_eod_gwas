# Sweep evidence at the six GWAS peaks: analysis design

This file describes the choices the evidence scripts implement (Supplementary Fig. S7, Table S5; ESM "Tests for
Selection at GWAS Peaks"). Code comments cite its section numbers. `../README.md` says how to run it.

## 1. Approach

Each sweep statistic is reduced by one fixed operator over a search region around each GWAS index SNP, and the
result is compared with several thousand non-peak regions reduced by the identical operator. The decision rule
is then applied to those same non-peak regions to measure how often it fires by chance. A peak window's maximum
is never compared with a distribution of single SNPs: a window holding hundreds of scored SNPs almost always
contains one above the single-SNP 99th percentile.

## 2. Units of analysis

Peaks are the extended GWAS peaks in `input_data/10_Introgression/peaks.bed` (SNPs with −log10 p > 3 merged
within 75 kb). Index SNPs are the GWAS index SNPs in column 6 of that file.

| Peak | Gene | Chromosome | Interval | Index SNP | Lineage(s) in which the BP allele arose |
|---|---|---|---|---|---|
| 1 | c6orf132 | chr6 | 29,448,847–29,509,385 | 29,503,232 | BP1 |
| 2 | itih3 | chr8 | 15,771,736–15,805,670 | 15,786,956 | BP2, BP3 |
| 3 | cebpg/cebpa | chr13 | 24,628,382–24,929,106 | 24,780,908 | BP2, BP3 |
| 4 | cntnap5b | chr16 | 10,728,257–11,082,159 | 10,916,883 | BP1 |
| 5 | sptbn4b | chr17 | 11,494,062–11,781,225 | 11,692,887 | BP1 |
| 6 | naa40 | chr24 | 8,331,219–8,730,509 | 8,596,847 | BP2 |

Lineages and their matched references (`input_data/11_Drift_vs_Selection/sweep_detection/split_up/`): BP1 = Bambomo
(n = 42) with TP1 = Ipassa + Bengue (n = 35); BP2 = Bikagala + Bavavela (n = 48) and BP3 = Biroundou (n = 16), both
with TP2 = Douvalou + Douengi (n = 31).

**Six tests, one per peak.** A peak's statistic is the maximum over its eligible lineages, and the null applies
the same maximum over the same lineages, so the multiple-lineage advantage is built into the null rather than
corrected for afterwards. BP2 and BP3 share BP alleles at Peaks 2 and 3, so they are not independent tests;
treating each peak as one test avoids claiming independence that does not exist. Per-lineage results are
reported alongside the six tests.

**Controls, reported but not part of the rule.** The matched TP population at each peak; the off-target BP
lineages; and BP1 at Peak 2, where Peak 2's BP allele did not arise, as a reminder that a high percentile alone
does not establish that a sweep drove the GWAS signal.

## 3. Applicability map

Step 1 (`01_index_snp_facts.R`) records, for every peak and lineage, the frequency of the BP-associated allele and
whether it is derived or ancestral. This determines which statistics are applicable:

- **iHS and nSL** compare haplotypes carrying the derived allele with those carrying the ancestral allele within
  one population and lose power as the allele approaches fixation. They are applicable only where the BP allele
  frequency is ≤ 0.85 in that lineage.
- **Rsb, XP-EHH, SweepFinder2 CLR, π ratio and PBS** do not require intermediate frequency and apply everywhere.

Where a statistic is not applicable it is scored `NA` and excluded from both the peak and the null, never scored
as a failure.

## 4. Statistics: one primary per evidence class

Within each class the available statistics are near-redundant (Rsb, XP-EHH, iHS and inES are built from the same
per-site integrated-EHH quantities), so one primary statistic per class keeps the composite rule in §6
interpretable.

| Class | Primary statistic | Window operator | Rationale |
|---|---|---|---|
| H — haplotype structure | Rsb (BP versus matched TP) | mean over the best 20-kb window | standardized, so it is not driven by the genome-wide BP/TP haplotype-length offset |
| S — site-frequency spectrum | SweepFinder2 CLR | maximum grid value in the region | detects sweeps after haplotypes have broken up; uses a genome-wide empirical background spectrum |
| D — diversity | π(BP)/π(TP) | minimum over the best 20-kb window | a sweep removes diversity in the swept lineage specifically |

**Secondary statistics**, reported in the evidence table but not part of the decision rule: XP-EHH; iHS
(`rehh::ihh2ihs`, freqbin 0.025); nSL (`scikit-allel`, standardized by allele-frequency bin), which is less
sensitive than iHS to variation in recombination rate; SweepFinder2 on the folded spectrum, as a check on
polarization error.

## 5. Search region, window operator, and null

**Primary configuration.** Search region: GWAS index SNP ± 100 kb. Window: the best (highest for Rsb and CLR,
lowest for the π ratio) 20-kb sub-window, stepping 5 kb, requiring at least 20 scored SNPs. For the CLR the
operator is instead the maximum grid value in the region, on the 10-kb grid.

The index SNP comes from the GWAS and is not chosen by any statistic computed here, so centring the search on it
does not use the sweep evidence twice.

**Operator identity.** The null applies the identical operator: same region size, window size, step, minimum SNP
count and maximum over lineages. The code asserts this (`assert_same_operator()` in `region_stats_lib.R`). The
alternatives in §7 are reported as separate rows, each with its own matched null, and never combined into the
headline by taking the best.

**Null sets.** Pseudo-peaks are positions every 200 kb on all 25 chromosomes, excluding every peak ± 250 kb.

- **(a) Unmatched pseudo-peaks.** How many regions of this size anywhere in the genome look as sweep-like as the
  peaks do.
- **(b) Matched pseudo-peaks, the headline null.** Matched by decile on SNP density and on callable fraction
  (pixy's per-window count of sites with data), relaxing the matching range until at least 1,000 regions are
  available. These covariates affect every statistic through the amount of information in a window, and neither
  is a consequence of a sweep.

**Not matched on differentiation.** A completed sweep creates high differentiation, so a null matched on
differentiation is built partly from other sweeps and from low-recombination regions where differentiation
accumulates by drift. Differentiation-matched results are therefore a sensitivity row (§7), where they bound how
much of the signal could be explained by the ascertainment of the peaks.

## 6. Decision rule and calibration

For each peak, each primary statistic is converted to a percentile against the matched null (b).

- **Supported:** at least two of the three primary statistics at or above the 95th percentile.
- **Strong:** supported, and at least one primary statistic at or above the 99th percentile.
- **Composite score:** the second-highest of the three class percentiles, the value the "two of three" rule turns
  on. The empirical p for a peak is the fraction of matched pseudo-peaks whose composite score is at least as
  high.

**Calibration.** The same rule is applied to every pseudo-peak in both null sets, giving the fraction of
non-peak regions called "supported" and "strong". The Spearman correlation matrix among all statistics across
the pseudo-peaks shows how much independent information the statistics carry.

**Lineage specificity is a label, not a gate.** Each peak is labelled *BP-specific* (matched TP below the 95th
percentile), *shared* (TP also at or above the 95th), or *unsupported*. The TP population may carry its own
sweep, and both populations rise together in regions of low recombination, so requiring the TP to be low would
discard real signals and, at n = 31 or 35, would mostly measure power.

## 7. Sensitivity rows (reported in full, never maxed into the headline)

Each is recomputed end to end, peak and null together: search region ± 50 kb; window size 40 kb; the extended
peak interval in place of a fixed search region and window; the unmatched null (a); a null also matched on BP–TP
differentiation; a null stratified by a recombination proxy (the mean inES of the lineages *not* in the contrast,
never the matched TP, which may share the signal under test); SweepFinder2 on the folded spectrum; XP-nSL in place
of Rsb; and every haplotype statistic recomputed with all lineages down-sampled to n = 16, the size of BP3, with a
fixed seed.

The down-sampled row has a control row. The only down-sampled Rsb is the one `10_genome_stats.R` recomputes from
the saved `scan_hh` output, standardized across all 25 chromosomes at once, whereas the headline Rsb is
standardized within each chromosome by `ines2rsb`. Comparing the down-sampled row directly with the headline would
vary the sample size and the standardization together. The control row (`_rsbgw`) is the full sample with the
genome-wide standardization, so `_rsbgw` against `_n16` isolates the sample size.

Down-sampling matters because the phasing (SHAPEIT4, no genetic map, no reference panel) has an unmeasured
switch-error rate that varies with sample size and heterozygosity, and switch errors shorten inferred haplotypes.

## 8. Interpretation

A peak that passes is consistent with a sweep on or near the BP-associated allele; it does not establish that the
sweep acted on the causal variant. A peak that fails has not been shown to be free of selection: a sweep that is
old, soft, or from standing variation, or one in a region of high recombination, may leave no detectable trace at
these sample sizes.

## 9. Outputs

In `output_data/11_Drift_vs_Selection/sweep_detection/evidence/`:

| File | Contents |
|---|---|
| `index_snp_facts.tsv`, `applicability_map.tsv` | BP-allele frequency and derived/ancestral state per peak per lineage |
| `evidence_matrix.tsv` | peak × statistic percentiles, composite score, empirical p, call, lineage-specificity label |
| `calibration.tsv` | false-positive rate of "supported" and "strong" in both null sets; per-statistic percentile uniformity |
| `statistic_correlations.tsv` | Spearman correlations among statistics across pseudo-peaks |
| `evidence_matrix_<row>.tsv`, `calibration_<row>.tsv` | each sensitivity row of §7 |
