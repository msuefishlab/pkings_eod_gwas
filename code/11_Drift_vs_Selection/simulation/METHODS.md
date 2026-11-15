# METHODS -- Peaks 4/5 simulation null (stages 00-07)

Exploratory. This file describes the coalescent null and the forward power grid run by `run_all.sh`. The
LD-decay (23*), empirical clustering (24*) and temporal ancestry (26*) scripts are described in their own
docstrings and in `README.md`.

## 1. Question

Peaks 4 (cntnap5b, chr16:10916883) and 5 (sptbn4b, chr17:11692887) are BP1/Bambomo peaks at which the
drift test (`demographic_modeling/fixed_differences/`) is uninformative by construction: its statistic E is
the expected number of genome-wide sites at least as BP-biased as the index SNP, and at Delta = 0.8084 and
0.7836 that lands in the top tier (E = hundreds to >1000) simply because the alleles are below fixation.
This analysis replaces that statistic with (a) a demography-informed coalescent null WITH LINKAGE, and (b) a
forward-simulation power arm with explicit dominance.

## 2. Fitted history carried forward

Pair A1w (whole-population Bambomo vs TP1 = Apassa + Bengue), from the drift test's
`output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences{,_agree}/fits/best_fits.tsv` and
`tail/models.tsv`.

The six-species-agreement polarization run (`fixed_differences_agree`, `RUN=agree`, the default) has seven
admissible A1w models (M0, M0f, M1f, M2, M2f, M3, M3f), which are used as they are. The rest of this section
describes the first-pass run (`fixed_differences`, `RUN=first`).

Four two-population models were fitted: M0 (split, no migration), M1 (split + symmetric migration),
M2 (split + asymmetric migration), M3 (isolation-with-migration, exponential growth). No folded variants
(M0f..M3f) exist for A1w in that run.

**All four are `converged = True` and `adequate = True` but `admissible = False`**, because the admissibility
rule in `09_tail.py` requires `p_misid <= 0.05` for unfolded analyses and A1w's fitted
ancestral-misidentification is 0.0785-0.0859 in every model. The drift test therefore reports A1w's tier as
"model-dependent (no admissible model): report O". Since the admissible set is empty, and since the drift
test's E_robust takes the maximum over admissible models, `RUN=first` **carries forward all four
converged-and-adequate models** and reports the maximum over them. That is the drift-generous choice, and it
departs from using only admissible models.

Convergence is judged in Godambe-adjusted log-likelihood units, as in the drift test
(`gap_adj = LL_gap_top3 / c_adj`: M0 0.000, M1 0.097, M2 0.015, M3 0.947; all < 1).

Unit conversion (recorded in `01_fitted_history_physical.tsv`):
`N_anc = theta / (4 mu L_mask)` with `L_mask = 532,739,484`, `L_scan = 618,475,302` (A1w.pkl),
`GENOME = 888,337,685`; `mu in {4e-9, 1e-8}` as in `09_tail.py` (`MU = (4e-9, 1e-8)`);
sizes `N_i = nu_i * N_anc`; time `T_gen = T * 2 N_anc`; migration `m_per_gen = m_dadi / (2 N_anc)`
with the dadi naming convention (`m_BP_from_TP` = forwards-time fraction of BP made up of TP migrants),
mapped to `msprime.Demography.set_migration_rate(source=BP, dest=TP, rate=m_BP_from_TP/(2 N_anc))`.

**Recombination rate: there is no recombination map for this species.** The sweep analysis's recombination
null (`code/11_Drift_vs_Selection/sweep_detection/evidence/06_recomb_null.sb`) is an empirical proxy-matched
null and assumes no per-bp rate; the drift test uses none. A rate therefore has to be chosen here and is an
assumption: baseline `r = 1e-8` per bp per generation (~1 cM/Mb, the order of magnitude reported for
teleosts), with the two mutation rates giving two linkage arms (r/mu = 2.5 at mu = 4e-9, r/mu = 1.0 at
mu = 1e-8). Sequence is simulated as the callable genome `L_scan` laid out contiguously, which mildly
OVERSTATES linkage among callable sites relative to their true spread over 888 Mb; overstated linkage
increases clustering under drift and is therefore the drift-generous direction.

## 3. Observed targets (from `index_units.tsv`, reproduced from the VCF)

A1w sample sizes: Bambomo 67 diploids (AN = 134), TP1 47 diploids (AN = 94).

| Peak | site | f_BP | f_TP | Delta | het_BP | het_TP | pooled f |
|---|---|---|---|---|---|---|---|
| 1 | chr6:29503232 | 134/134 | 1/94 | 0.9894 | 0 | 1 | 0.5921 |
| 4 | chr16:10916883 | 124/134 | 11/94 | 0.8084 | 6 | 9 | 0.5921 |
| 5 | chr17:11692887 | 105/134 | 0/94 | 0.7836 | 23 | 0 | 0.4605 |

Genotype counts are re-extracted from the all-sites VCF for the 114 A1w fish by
`00_extract_index_genotypes.sh` and must reproduce every `index_units.tsv` allele count exactly
(`02_observed_targets.tsv`).

## 4. Statistic set and decision rule

For each model in {M0, M1, M2, M3} x each mu in {4e-9, 1e-8}, simulate replicate callable genomes under
`msprime` with the fitted history, sampling 134 Bambomo + 94 TP1 chromosomes, and record for every
segregating site the joint tuple (d_BP, d_TP, het_BP, het_TP). Statistics, per replicate genome:

- **S1 (marginal)** `N_delta(x)` = number of sites with Delta >= x, for x = 0.8084 (Peak 4) and 0.7836
  (Peak 5), polarised in the BP-biased direction only (one-sided, as in the drift test).
- **S2 (joint)** `N_joint` = number of sites with Delta >= x AND the observed heterozygote configuration
  (het_BP, het_TP) exactly.
- **S3 (clustering, the linkage-specific statistic)** the maximum number of Delta >= x sites inside any
  500-kb window of a simulated genome, compared with the observed count inside the Peak 4 and Peak 5
  regions of `peaks.bed`.

Reported as mean per genome with the 2.5-97.5% range over replicates, and the **maximum over the models**
as the drift-generous summary, in direct analogy with the drift test's `E_robust`.

**Decision rule.** Tiers are read on the same ladder as the drift test's `09_tail.py` `tier_of()`:
E <= 0.005 "drift implausible (all units)"; <= 0.05 "drift implausible (this unit)"; <= 1 "unusual under
drift"; <= 20 "drift makes a few"; > 20 "drift readily makes this". The tier is read from the
drift-generous (maximum-over-model) value of each statistic, separately for Peak 4 and Peak 5 and separately
for S1, S2 and S3.

**Calibration precondition.** The joint statistic S2 is only interpretable if the neutral simulations
reproduce the observed heterozygote baseline. The precondition is tested as: do the observed heterozygote
counts at Peaks 1, 4 and 5 fall inside the central 95% of the simulated het distribution CONDITIONAL on the
observed per-population allele counts? If they do not, S2 is declared **uncalibrated** and only S1 and S3
are reported. If the observed hets are at Hardy-Weinberg expectation within each population, the observed
heterozygote contrast between the populations is a consequence of the allele frequencies rather than
independent information, and S2 then cannot be more extreme than S1 by more than the HWE factor.

## 5. Power arm (SLiM)

Forward simulations under the same fitted history (the single most drift-generous / highest-LL model,
rescaled), with a single selected site. Grid:
- h in {0, 0.5} for Peak 4 (recessive / additive coding) and h = 1 for Peak 5 (dominant coding).
- s in {0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1}.
- onset in {0.1, 0.25, 0.5, 1.0} x (generations since the split), and starting frequency
  {de novo 1/2N, standing 0.01, standing 0.05}.
A (h, s, onset, f0) cell is **compatible** if the BP-population frequency of the selected allele at
sampling falls in the observed binomial 95% interval for d_BP/AN_BP AND the TP frequency falls in the
observed interval for d_TP/AN_TP. Compatibility fraction over replicates is reported; cells with
compatibility fraction >= 0.05 define the compatible region.

**Power floor.** The grid's smallest s and most recent onset define what this design could NOT have
detected, and this is reported quantitatively (the s below which the compatible region is
indistinguishable from the neutral replicates).

The agreement-run grid (`05d_gen_slim_tasks_agree.py`) uses the seven admissible histories, h in
{0, 0.5, 1}, s down to 0.0002 and onsets of 25-5,000 generations, with neutral controls included.

## 6. Not addressed

Standing variation recruited at the moment of the split (f0 = ancestral polymorphism at arbitrary
frequency), polygenic / many-route architectures, selection acting on a linked site rather than the index
SNP, and any history outside the fitted two-population models.

## 7. Additional statistics

1. **Window-matched clustering.** S3 is computed on the simulation's own 2-Mb chunks
   (`08_clustering_2Mb_windows.tsv`) and, from the recorded positions of qualifying sites
   (`04b_worker_msprime_positions.py`, `04c_worker_msprime_agree.py`), in windows matched to the peak spans
   (353,902 and 287,163 bp; `15_clustering_matched_window.tsv`). The matched-window version is the headline.
2. **Neutral forward controls.** Cells with s = 0 at every (model, f0, onset) and 500 replicates
   (`05c_gen_slim_tasks_controls.py`; `13_neutral_forward_controls.tsv`), plus a finer s grid
   {0.0002 ... 0.003} at onsets 1000-5000 to locate the power floor (`12_power_grid_combined.tsv`,
   `14_power_floor.tsv`).
3. **Empirical window baseline.** S3 compares a GWAS-ascertained region against a simulated null, and the
   regions were defined by association with the phenotype, which is itself correlated with BP-TP
   differentiation. To bound that ascertainment, the empirical genome-wide distribution of 354-kb and
   287-kb window counts is computed from the same site table (`03_window_scan.py`; `16_*`, `17_*`). Where
   the simulated and empirical window nulls disagree, the empirical one is the more conservative and is
   reported as such.
