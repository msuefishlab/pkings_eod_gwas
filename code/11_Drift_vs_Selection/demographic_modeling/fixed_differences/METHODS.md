# Drift test at the GWAS index SNPs: methods as implemented

This file describes what the scripts in this directory compute. The paper's account is the ESM section "Drift Test
at GWAS Peaks" and the Table S6 legend. Section numbers here are the ones the script comments cite
("section 4", "section 9.3", ...).

## 1. Question

For each GWAS index SNP in each BP population: how many SNPs at least as differentiated as the index SNP would
neutral drift produce anywhere in the genome, under a two-population history fitted to the genome itself (E), and
how many are observed outside the GWAS peaks (O)?

## 2. Reading E

- E is an expected count, so by Markov's inequality P(drift produces at least one such SNP) ≤ E, whatever the
  linkage among sites.
- The fitted history includes old polymorphisms that sort differently in the two daughter populations, so E is
  the null for "drift sorted standing variation", not only for new mutations.
- A large E means differentiation of this size is compatible with drift; it does not show that drift raised the
  index allele.
- **K_needed** = L_scan / E_genome (best model): how many BP-capable sites the genome would need for drift to be
  expected to produce one this differentiated. `13_report_tables.py` also gives 1 / E_robust, the share of all
  callable sites that would each have to be BP-capable.

## 3. Units, samples and data

**Pairs** (`pairs.tsv`, written by `00_prepare_inputs.py`; BP is population 1):

| Pair | BP | TP reference | Own index SNPs |
|---|---|---|---|
| A1w | Bambomo (BAM) | TP1 = Apassa + Bengue | Peaks 1, 4, 5 |
| A2w | Bikagala (BIK) | TP2 = Douvalou + Douengi | Peaks 2, 3, 6 |
| A3w | Bavavela (BAVA) | TP2 | Peaks 2, 3, 6 |
| A4w | Biroundou (BIR) | TP2 | Peaks 2, 3 |
| C1w | Apassa | Bengue | TP–TP control |
| C2w | Douvalou | Douengi | TP–TP control |

The 11 own units are the BP pair × own index SNP combinations. Every pair is also evaluated at the other index
SNPs (off-target controls). Index SNPs are column 6 of `input_data/10_Introgression/peaks.bed`.

**Groups** (`groups.tsv`, `group_sizes.tsv`): `w*` groups are whole populations from
`input_data/05_PopGen/popsfile.txt`, mixed-phenotype fish included (wBAM 67, wAPA 30, wBEN 17, wBIK 30, wBAVA 23,
wBIR 16, wDOV 16, wDOG 15); these are the reported pairs. `p*` groups are the phenotype-pure GWAS lists
(`input_data/11_Drift_vs_Selection/sweep_detection/split_up/`), used for sensitivity pairs and for the index-SNP
check of section 9.1. One `OUT_*` group per outer-clade outgroup species (section 3, Ancestral allele).

**Data.** The unphased all-sites VCF (`${allsites_vcf}`; SNPs and invariant sites), read once per 10-Mb chunk
(`chunks.tsv`) by `01_site_table.sb`:
- genotypes with no reads (FMT/DP < 1) are set to missing before AN/AC are counted (GATK writes 0/0 for
  uncovered samples at invariant sites; pixy does not count them);
- per-group AN, AC and ExcHet from `bcftools +fill-tags`;
- **paralog filter:** SNPs with ExcHet < 0.001 in any whole population are dropped from the data and from L;
- **callable site for a pair:** AN ≥ the projection size in both populations; invariant sites must also be
  outgroup-callable. Only biallelic SNVs on chr1–25 are used.

**Ancestral allele.** The six outer-clade species are separate outgroup groups: *P. szaboii* (OUT_PSZA, 11 fish),
*P. magnostipes* (OUT_PMAG, 20), *P. parsegnis* (OUT_PPAR, 7), *P. hopkinsii* (OUT_PHOP, 5; PSZA_455 is
*P. hopkinsii* and is counted here), *P. gabonensis* (OUT_PGAB, 4) and *P. vadamans* (OUT_PVAD, 2). A SNP is
polarized only if at least 2 species have data (AN ≥ 2 each), they carry at least 8 alleles in total, and every
species with data is fixed for the same allele; otherwise it is unpolarized. An invariant site is
outgroup-callable under the same data rule. The *P. curvifrons* group is not used, because it is closer to
*P. kingsleyae* and shares excess derived alleles with it.

**Two analyses per pair.**
- *Unfolded:* polarized SNPs, models with an ancestral-misidentification parameter, derived-allele Δ.
- *Folded (label-free):* every SNP, polarized or not, ALT-oriented, fitted folded with models M0f–M3f (no
  misidentification parameter). Its tail counts SNPs at least as differentiated in either direction, so its E can
  only exceed the directional E.
- *Headline per unit:* the unfolded analysis if the pair has an admissible unfolded model and the index SNP is
  polarized; otherwise the folded analysis. The Peak 2 and Peak 6 index SNPs are unpolarized (the outgroups
  disagree), so their units use the folded analysis. Both analyses are reported.

**Masks used for fitting** (`masks/exons.bed`, `masks/exclude.bed`):
- exons excluded (all transcript types, RefSeq GFF); exons ±5 kb is a sensitivity configuration
  (`06_gather.py --mask_c3 2`);
- the six extended GWAS peaks ±250 kb and the chr4 inversion-like polymorphism (chr4:40.88–42.77 Mb) excluded;
- **call deserts** (`05_deserts.py`): per pair, a 100-kb bin is excluded when its callable SNPs polymorphic in the
  pair, per callable site, fall below 10% of the genome-wide median for that pair. Desert bins are removed from
  L_mask and L_scan alike.

**Projection** (`04_projections.py`): hypergeometric. Per population, n = the 10th percentile of AN over masked,
polarized SNPs, rounded down to an even number and capped at 60 haplotypes (`--cap 40` lowers the cap if model
evaluation is too slow).

## 4. Statistics

Per unit (BP population, TP reference, index SNP); i and j are derived-allele counts in BP and TP.

- **Δ_idx:** derived-allele frequency in BP minus that in TP at the index SNP, from raw counts (folded analysis:
  the ALT-allele difference, used as |Δ|).
- **Tail region R:** projected cells with i/n_BP − j/n_TP ≥ Δ_idx (folded: |i/n_BP − j/n_TP| ≥ |Δ_idx|). A
  tolerance of 10⁻⁴ keeps the index SNP's own cell in R despite rounding of Δ.
- **O:** observed SNPs in R; O_mask over the fitting mask, O_genome over all callable sites outside the unit's own
  peak ±250 kb.
- **E_mask:** the fitted model's θ-scaled expected count in R. **E_genome** = E_mask × L_scan / L_mask, where
  L_scan is every callable site on chr1–25 outside call deserts.
- **E_hetNe** (local effective size): per 1-Mb block w, f_w = (π / d_OUT)_w over the whole populations not in the
  pair, divided by the same ratio genome-wide (d_OUT, divergence to the outgroup); blocks whose outgroup
  divergence is below 25% of the median get f = 1. The per-site tail probability q is recomputed at local size
  N × f (dadi units T/f, m·f, θ·f) on a grid of f, interpolated with f clipped to 0.2–1.5, and
  E_hetNe = Σ_w L_w q(f_w), with L_w the block's callable sites over the whole scan (`06b_block_scan.py`). The same
  sum on f in 0.05–2.0 is reported as information. **f_peak** is the index SNP's own block value; **f_star** the
  largest f at which a 500-kb peak region alone would be expected to hold one such site.
- **τ_BP, τ_TP:** each branch's total drift since the split, ∫dt / 2N(t) from the fitted model; F = 1 − e^(−τ);
  T/2N_anc.
- **Ladder:** O and E over the fitting mask at Δ ≥ 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0.
- **Asymmetry** (unfolded models only): observed and expected BP-derived (Δ ≥ 0.95) and TP-derived (Δ ≤ −0.95)
  near-fixed differences.
- **Parallel count** (`11_clusters.py`; Peaks 2, 3, 6): SNPs outside every extended peak where one allele is at
  least as frequent as at the index SNP in each BP population of the lineage and at most as frequent in TP2.

## 5. Models and fitting

- dadi 2.4.3 (`images/dadi_tools.sif`), grid pts = [n+10, n+20, n+30] with extrapolation. BP is population 1, so
  in M3 the founder fraction s belongs to BP.
  - **M0** split, no migration; **M1** split, symmetric migration; **M2** split, asymmetric migration;
    **M3** isolation with migration with a founder fraction s and exponential growth.
  - Unfolded models carry dadi's ancestral-misidentification parameter; M0f–M3f are the same demographies fitted
    to the folded label-free spectrum.
- **Starts** (`run_fits.sh`, `07_fit_models.py`, BOBYQA via `dadi.Inference.opt`): round 1, 20 random perturbed
  starts per pair × model; round 2, 20 starts perturbed (fold 0.5) from the best start so far; round 3, 20 starts
  (fold 0.25) from the best start so far.
- **Nested-start guard:** M0 ⊂ M1 ⊂ M2 within each analysis. Where a model's best log-likelihood lies below a
  sub-model's, it gets 20 more starts from the sub-model's best fit embedded in it (M0 → M1: m = 0; M1 → M2: both
  rates = m), perturbed by fold 0.25 (`run_fits.sh nested 1` after round 2, `nested 2` after round 3). A model
  still below a sub-model at collection is not converged (`nested_ok`). M3 does not nest M2 and is not checked.
- **Convergence** (`08_collect_fits.py --godambe`): the best three starts, pooled over all rounds, lie within one
  Godambe-adjusted log-likelihood unit, (LL₁ − LL₃)/c ≤ 1, where c is the mean eigenvalue of H⁻¹J at the best start
  (the factor by which linkage inflates the composite log-likelihood). A non-positive c means not converged. The
  raw rule (best three within 1 unit, parameters within 5%) is reported as `converged_strict`.
- **Uncertainty:** Godambe (GIM) from 100 bootstrap spectra over 1-Mb blocks (`--block_mb 5` for 5-Mb blocks).
  Godambe-adjusted LRTs (M0 ⊂ M1 ⊂ M2) are reported for information; no AIC (composite likelihood).
- **Ridge rule:** for a converged model, E is the largest E among its best three starts.
- **Admissible model:** converged, misidentification rate ≤ 0.05 (unfolded), and passing the adequacy gate: at the
  most extreme ladder threshold with E ≥ 10, 1/3 ≤ O/E ≤ 3.
- **Best model:** the admissible model with the highest log-likelihood, per analysis. 100 block-bootstrap refits
  under the best model and under M3 (M3f) give 95% intervals for E, T/2N and τ; 20 Poisson replicates of the best
  model test parameter and E recovery (section 9.3).

## 6. Decision rule

### 6.1 E_robust

The largest of: the upper 97.5% bootstrap bound of E_genome under the best model; the point E_genome under every
admissible model; E_hetNe under the best model; and E_grid2, E under the best model on a doubled grid
(pts → 2 pts, θ refitted).

### 6.2 Tiers

| E_robust | Reading |
|---|---|
| ≤ 0.005 | drift implausible across all units |
| 0.005–0.05 | drift implausible for this unit |
| 0.05–1 | unusual under drift |
| 1–20 | drift makes a few such sites |
| > 20 | drift readily makes this |

### 6.3 Overrides (as applied by `09_tail.py`)

- **Model fails at the extremes:** if E_robust ≤ 1 but O_genome > max(5, 3 × E_genome), the tier is read from
  O_genome.
- **No admissible model for an analysis:** no E is assigned; O_genome is reported.
- **Model-dependent:** a unit whose tier differs between admissible models is tagged "model-dependent".
- **Numerically unreliable:** a best model flagged by section 9.3 (corner cells move > 5% on the doubled grid) is
  tagged.

### 6.4 Drift exposure

Reported per pair next to the tiers (`tail/drift_exposure.tsv`): whether τ_BP > τ_TP under the best model, and the
model-relative asymmetry of near-fixed differences (section 4).

### 6.5 Controls

Off-target index SNPs in every pair and the two TP–TP pairs, which should show O/E near 1 at the corners.

## 7. Sensitivity configurations

Supported by the scripts through `--pairs_role`/`--suffix` (03, 06, 07, 08, 09), `--mask_c3` (06) and
`--block_mb` (08): phenotype-pure groups (p*), single TP references (APA, BEN, DOV or DOG), BIK + BAVA pooled
(B2), exons ±5 kb mask, and 5-Mb bootstrap blocks.

## 8. Near-fixed differences and their clusters

`11_clusters.py` counts BP-derived and TP-derived near-fixed SNPs (≥ 0.95 in one population and ≤ 0.05 in the
other; also 0.9 / 0.1) outside the extended peaks, joins those within 75 kb into clusters, and writes the cluster
BEDs, `fixed_diff_cluster_counts.tsv` and `parallel_count.tsv`. `phaseB_gate.tsv` ("Phase B" is the cluster null
for the sweep statistics in `11_Drift_vs_Selection/folded_clr/`) gives, for the lineage pairs pBAM/pTP1, pBP2/pTP2
and pBIR/pTP2, the threshold whose clusters form that null: 0.95 / 0.05 if the lineage has at least 200 non-peak
clusters there, else 0.9 / 0.1 if it has 200 there, else none.

## 9. Checks

### 9.1 Site table (`02_check_site_table.py`; downstream jobs run only if it passes)
- whole-population π per population and chromosome, and dxy for six check pairs, within 0.5% of pixy
  (`pixy_reference_totals.tsv`);
- SNP record count reported against the phased SNP set (not gated);
- index-SNP ALT frequencies in the phenotype-pure lineage groups within 0.05 of
  `sweep_detection/evidence/index_snp_facts.tsv`;
- the BP allele not called ancestral by the outgroups (unpolarized index SNPs are reported as information).

### 9.2 Spectra (`10_checks.py`)
Block spectra sum to the genome-wide spectrum; spectrum total plus invariant sites equals L_mask; the folded fit
target equals the fold of the label-free spectrum, and the folded polarized spectrum lies cell by cell within the
folded label-free one.

### 9.3 Models (`10_checks.py`, `09_tail.py`)
Model-implied π, dxy, Fst and Tajima's D against the data (the sign of Tajima's D must match on the BP side for
BAM and BIK); corner cells (|Δ| ≥ 0.95) must change by less than 5% when the grid doubles, otherwise E is tagged
"numerically unreliable"; Poisson parameter and E recovery from 20 replicates of each best model.

## 10. Outputs

In `${FD_OUT}` (`output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree/`):
`site_table/`, `checks/`, `sfs/` (spectra, projections, deserts, block scans), `fits/` (`best_fits.tsv`,
`lrt.tsv`, per-start files), `tail/` (`models`, `ladder`, `unit_results`, `drift_exposure`, `recovery`, `curve`),
the cluster files of section 8, the tail figure (`12_tail_figure.R`) and the report tables
(`13_report_tables.py`), including one sentence per own unit: since the split drift has moved a typical 50% allele
by about ±√(F·0.25); the elapsed time is T/2N_anc and, in the BP population's own drift units, τ_BP, given as a share
(τ_BP / 2) of the time a new neutral mutation needs to drift to fixation; under the best history drift would create
E such sites (bootstrap 95% interval; E_robust); O are observed outside the GWAS peaks; tier.
