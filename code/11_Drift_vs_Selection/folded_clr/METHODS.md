# METHODS: Peak 6 folded CLR and the fixed-difference null

Exploratory. Part A gives the "CLR on the folded spectrum" sensitivity row of Table S5 and the six-test
Benjamini–Hochberg correction quoted in the Fig. S7 legend. Part B scores Peak 6 against a null built from
drift-made fixed differences.

---

## Part A: folded-spectrum SweepFinder2 CLR

### A1. Question
Peak 6 (naa40, chr24:8,596,847, BP2 = Bikagala + Bavavela) has its composite score set by the SweepFinder2 CLR
(the second-highest of its three primary percentiles). The CLR on an unfolded spectrum depends on the
ancestral-allele call, and Peak 6's ancestral state is uncertain (*P. vadamans* carries the BP allele). The
folded CLR removes that dependence. Question: is Peak 6's status robust to polarization?

### A2. Input construction
`sweep_detection/sweepfinder2/vcf_to_sweepfinder2.py --polymorphic-only --folded`. It emits SweepFinder2
column 4 = 1 and x = min(d, AN − d), where d is the AA-polarized derived count.

Site set: identical to the unfolded run, i.e. still restricted to sites carrying a usable INFO/AA call and still
`--polymorphic-only`. Holding the site set fixed makes the folded-versus-unfolded difference attributable to
polarization alone. Consequence: this does not recover sites that the unfolded run dropped for lack of an AA
call, so it is a test of *mis*-polarization, not of AA coverage.

### A3. Correctness precondition (run before scaling up)
On chr1 for lineage BP2, from the same polarized VCF, require all of:
  P1. folded x == min(AC, AN − AC) at 100% of emitted sites (the folded count is invariant to the AA call);
  P2. folded x == min(d, AN − d), where d is the unfolded run's derived count, at 100% of sites present in both;
  P3. the emitted position set is identical between the folded and unfolded chr1 files;
  P4. the number of sites with AA == ALT is > 0 (otherwise folding is vacuous and the check uninformative);
  P5. column 4 is 1 on every folded row and 0 on every unfolded row.
If any of P1–P5 fails, the folded scan is not run. Output: `audit/folded_conversion_check.tsv`.

### A4. Scan
SweepFinder2 from `images/sweepfinder2.sif`, `-lg 10000` (10-kb grid, as in the unfolded scan), genome-wide
chr1–chr25, for all five lineages (BP1, TP1, BP2, TP2, BP3), because the recalibration in A5 scores all six peaks
and the null. The per-lineage empirical background spectrum is recomputed folded (`SweepFinder2 -f` over the
concatenated 25 chromosomes); the unfolded spectrum is not reused. The merge uses the column layout of
`sweep_detection/sweepfinder2/13_sf2_merge.sh`, written to `merged/<pop>.genome_wide_clr_folded.txt`, the file
that the `clr_folded` entry of `02_region_stats.R` reads.

BP2 is Bikagala + Bavavela pooled, as in the headline analysis; the two are not scanned separately because no
unfolded counterpart exists for a split BP2.

### A5. Recalibration with the headline operator
`sweep_detection/evidence/02_region_stats.R --stats rsb,clr,pi_ratio,clr_folded`, i.e.
  - search region = index SNP ± 100 kb;
  - CLR (folded and unfolded) reduced by `region_max`, min 3 scored grid points;
  - Rsb by best 20-kb sub-window mean, 5-kb step, min 20 sites, direction high;
  - π ratio by best 20-kb sub-window mean, min 4 bins, direction low;
  - null = the deterministic pseudo-peak grid (centres every 200 kb on chr1–chr25, minus every GWAS peak
    ± 250 kb), reduced by the same operator;
  - matched null = same SNP-density decile and callable-fraction decile, caliper widened one decile at a time
    until ≥ 1,000 pseudo-peaks qualify.
Re-deriving the unfolded `clr` column in the same run is a self-check: it must reproduce the headline
`evidence_matrix.tsv` percentiles for all six peaks (tolerance 1e-6 on `peak_value` and `pct_matched`), or
nothing downstream is reported.

### A6. Decision rule
The headline rule with `clr_folded` substituted for `clr` in the primary triple:
  - a peak is **supported** if ≥ 2 of its 3 primaries reach the 95th percentile of its matched null; **strong** if
    additionally ≥ 1 reaches the 99th;
  - composite = the second-highest of the three percentiles;
  - empirical p = fraction of matched pseudo-peaks whose own composite ≥ the peak's;
  - six-test correction over the six peaks: Bonferroni = 6p, Benjamini–Hochberg over the same six.
Verdict for Peak 6:
  - **ROBUST**: still supported and BH q ≤ 0.05;
  - **WEAKENED**: still supported but BH q > 0.05;
  - **NOT ROBUST**: no longer supported.
Reported for all six peaks, folded beside unfolded.

---

## Part B: the fixed-difference null

### B1. Precondition (≥ 200 drift-made clusters per lineage)
Read from `fixed_diff_cluster_counts.tsv` in
`output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree/` (columns pair / threshold /
direction / n_snps / n_clusters / n_snps_delta / an_rule). Pairs: A1w = BAM vs APA+BEN (BP1, Peak 1),
A2w = BIK vs DOV+DOG, A3w = BAVA vs DOV+DOG (both BP2, Peak 6). BP-derived, threshold 0.95. The counts are
checked against the line counts of the corresponding `fixed_diff_clusters.<pair>.bed`. A lineage with < 200
clusters fails the precondition for Part B at that lineage, and no substitute null is run.

### B2. Null construction
Null loci = the BP-derived near-fixed-difference clusters at Δ ≥ 0.95 for A2w and A3w, centred on BED column 5
(`centre_pos`, the most differentiated SNP of the cluster). The clusters are built outside every extended peak
± 250 kb by `demographic_modeling/fixed_differences/11_clusters.py`, the same buffer as the sweep null, so no
further exclusion is applied. A2w and A3w clusters are pooled into one BP2 null (matching the pooled BP2 lineage
that carries Peak 6's signal), and also reported separately.

### B3. Statistic and percentile
Each cluster centre is reduced by the operator of A5, for Rsb, CLR (unfolded), CLR (folded) and the π ratio, with
the same over-lineage aggregation as Peak 6 (lineage set BP2). Peak 6's percentile in this null = fraction of
cluster loci whose value is at or below Peak 6's, in the sweep-like direction (reversed for the π ratio).
Clusters whose region fails the operator's minimum-site rule return NA and are excluded from the denominator,
which is reported.

### B4. Decision rule
Peak 6 exceeds the drift-made null for a statistic if its percentile in that null is ≥ 95. Overall, Peak 6
survives Part B if ≥ 2 of the 3 primaries (Rsb, CLR, π ratio) reach the 95th percentile of the drift-made null,
with the CLR taken in whichever polarization Part A leaves standing; both polarizations are reported.

### B5. What this null is and is not
Every null locus is a place where drift drove a difference to near fixation in these same populations, so a
sweep statistic that cannot separate Peak 6 from them has not separated selection from drift. The clusters are
not matched to Peak 6 on SNP density, callable fraction or recombination, so a percentile in this null confounds
those with the drift history. Pooled and per-lineage percentiles are both reported.

---

## Scope
- BP2 is scanned as the pooled Bikagala + Bavavela population (A4).
- Folded inputs keep the unfolded run's AA-filtered site set (A2).
- `02_region_stats.R` is run with `--stats rsb,clr,pi_ratio,clr_folded` only; xpehh, ihs and nsl are secondary
  statistics and are not recomputed.
- Only the headline `matched` null mode is recomputed folded; the other sensitivity rows are not re-run folded.
