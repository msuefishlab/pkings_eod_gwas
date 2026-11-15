#!/usr/bin/env python3
"""
TFBS scan, Stage B: matched-null test of directional TFBS changes at the GWAS peaks.

Specification: code/15_TFBS_Scan/README.md.

For every peak unit (one SNP, or a cluster of SNPs <= F bp apart) the BP and TP haplotypes are compared motif by motif
(BP-only / TP-only calls, FIMO p < 1e-4). Each peak unit is matched to random non-coding pool units of the same
sequence composition (stratum), oriented by the rule in README.md, 10,000 times per locus and set. All statistics are
compared to that matched null with the two-sided empirical p of README.md.

Inputs (formats in README.md; $OUT = output_data/15_TFBS_Scan):
  --units      $OUT/units.tsv                (build_units.py)
  --hits       $OUT/hits.npz                 (parse_fimo_hits.py)
  --families   input_data/15_TFBS_Scan/jaspar2024_families.tsv
  --polarized  $OUT/peak_snps_polarized.tsv  (optional: cross-checks bp_mask, adds delta_af to the call tables)
Outputs in --out_dir ($OUT/results): summary.txt, locus_tests.tsv, family_tests.tsv, named_group_tests.tsv,
  balance.tsv, calibration.tsv, fallbacks.tsv, per_unit_calls_{locus}.tsv, null_fraction_{locus|pooled}.png

Usage:
  python3 code/15_TFBS_Scan/matched_null.py \\
      --units output_data/15_TFBS_Scan/units.tsv \\
      --hits output_data/15_TFBS_Scan/hits.npz \\
      --families input_data/15_TFBS_Scan/jaspar2024_families.tsv \\
      --polarized output_data/15_TFBS_Scan/peak_snps_polarized.tsv \\
      --out_dir output_data/15_TFBS_Scan/results

Choices that README.md leaves open (also printed in summary.txt):
  - f, events_per_snp and the sign-flip test count motif-level calls; family and named-group statistics count
    family-level (collapsed) calls. A family-level f is reported as a secondary statistic (`f_family`).
  - Null events_per_snp divides by the SNP count of the drawn pool units (a k>KMAX peak unit draws a KMAX chain).
  - k>=2 orientation: pool SNP j is matched to peak SNP j (j < pool k, so the first KMAX peak SNPs when peak k > KMAX).
    For C>G and T>A (self-complementary) both orientations reproduce the class, so bp_is_ref is copied.
  - Strata: GC quintile edges come from pool singletons for k = 1 and from pool chains for k >= 2; af edges from pool
    singletons. k >= 2 strata that hold < min_stratum chains drop gc, then gap_bin (not specified in README.md).
    If no level reaches min_stratum, the coarsest non-empty level is used and flagged `below_min`.
  - "Without replacement within a draw" holds within each group of peak units sharing a stratum (vectorised Floyd
    sampling); collisions across different strata of one draw are allowed.
  - Families/groups enter BH (or Bonferroni) only when testable: observed net != 0 or some null draw != 0.
  - The `Unknown` family is excluded from the recurrence count (it is not a family).
  - Power check: at each locus the family with the median null SD of nBP-nTP (extended set, standard rule) gets +20
    BP-only calls on top of the observed peak; reported as q. Also reported: the share of testable families that
    reach q < 0.05 the same way, and the power over 100 pseudo-peaks (null draws + 20).
  - SD guard aborts (exit code 2) only for the primary analysis (set extended, standard rule, every locus and pooled);
    other sets/rules only warn. calibration.tsv is written before the abort.
"""

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats

# --------------------------------------------------------------------------------------------------------------------
# Constants (README.md)
# --------------------------------------------------------------------------------------------------------------------
F = 32                       # flank, and max gap between consecutive SNPs of one unit
PRESENT_MLOG10P = 4.0        # presence call: FIMO p < 1e-4  <=>  -log10 p > 4
SCAN_MLOG10P = 3.0           # any hit:       FIMO p <= 1e-3 <=> -log10 p >= 3 (strong-change rule)
MLOG10P_TOL = 1e-6
SD_GUARD_FACTOR = 0.8
PLANT_N = 20                 # synthetic BP-only calls for the power check
LOCUS_ORDER = ["sptbn4b", "cntnap5", "itih3"]
SETS = ["extended", "core"]
RULES = ["standard", "strong"]
GAP_BINS = ["1-4", "5-12", "13-32"]
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
BASES = "ACGT"
COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
PRIMARY = ("extended", "standard")
# units.tsv column `af_ref` holds GEMMA `af`, which tracks the ALT allele frequency (checked at 8 SNPs).
# It is used only for the af-bin quintiles (peak and pool hold the same quantity); never call it the REF frequency.
AF_LABEL = "gemma_af (~ALT freq)"

# Named groups (statistic 4b): the BCL6 and C/EBP motifs. JASPAR files BCL6 inside a 124-motif zinc-finger family,
# so the family test cannot answer for it. BCL11A/BCL11B are unrelated motifs and are excluded; the PAR-bZIP members of
# JASPAR's "CEBP-related" family are left to the family test. The groups were fixed before the test was run.
NAMED_GROUPS = {
    "BCL6": ["MA0463.3", "MA0731.1"],  # BCL6, BCL6B
    "CEBP": ["MA0102.5", "MA0466.4", "MA0836.3", "MA0837.3", "MA0838.1", "MA1636.2", "MA0019.2"],  # C/EBP proper
}

UNIT_COLS = ["unit_idx", "source", "locus", "set", "chr", "k", "pos", "ref", "alt", "snp_off", "seq_start", "seq_len",
             "gap_bin", "gc", "cpg", "af_ref", "snp_id", "log_p", "bp_mask"]


def log(msg):
    print(f"[matched_null] {msg}", flush=True)


class CalibrationError(RuntimeError):
    """The matched null failed the SD guard (README.md statistic 6)."""


# --------------------------------------------------------------------------------------------------------------------
# Mutation classes
# --------------------------------------------------------------------------------------------------------------------
def directed_class(bp, tp):
    """BP->TP base change, strand-collapsed so the BP base is a pyrimidine (C or T)."""
    bp, tp = bp.upper(), tp.upper()
    if bp not in COMPLEMENT or tp not in COMPLEMENT or bp == tp:
        raise ValueError(f"invalid base change {bp}>{tp}")
    if bp in "AG":
        bp, tp = COMPLEMENT[bp], COMPLEMENT[tp]
    return f"{bp}>{tp}"


BASE_CODE = {b: i for i, b in enumerate(BASES)}
DC_TABLE = np.full((4, 4), -1, dtype=np.int64)          # DC_TABLE[bp_code, tp_code] -> class index
for _b in BASES:
    for _t in BASES:
        if _b != _t:
            DC_TABLE[BASE_CODE[_b], BASE_CODE[_t]] = CLASSES.index(directed_class(_b, _t))
# reverse of a class: the same SNP read in the other orientation (C>A <-> T>G, C>T <-> T>C, C>G and T>A fixed)
REV_CLASS = np.array([CLASSES.index(directed_class(c[2], c[0])) for c in CLASSES])
UND_CLASS = np.minimum(np.arange(6), REV_CLASS)          # undirected class id
SELF_COMP = REV_CLASS == np.arange(6)                    # C>G, T>A


def base_codes(arr):
    """Array of single-base strings -> int codes 0..3; raises on anything but ACGT."""
    arr = np.asarray(arr).astype(str)
    lut = np.full(256, -1, dtype=np.int64)
    for b, i in BASE_CODE.items():
        lut[ord(b)] = i
        lut[ord(b.lower())] = i
    bad = np.char.str_len(arr) != 1
    if bad.any():
        raise ValueError(f"non-single-base allele(s), e.g. {arr[bad][:3].tolist()}")
    codes = lut[np.frombuffer("".join(arr.tolist()).encode("ascii"), dtype=np.uint8)]
    if (codes < 0).any():
        raise ValueError(f"alleles must be ACGT, found {sorted(set(arr[codes < 0].tolist()))[:5]}")
    return codes


def orient_chain_masks(pool_ref, pool_alt, peak_bp, peak_tp, peak_bits):
    """Orientation (BP mask) of pool chains matched to one peak unit (README.md, k >= 2).

    pool_ref, pool_alt: (n_cand, kc) base codes. peak_bp, peak_tp, peak_bits: per peak SNP (position order); only the
    first kc are used. For SNP j: if the pool SNP's undirected class equals peak SNP j's, choose the orientation that
    reproduces peak SNP j's directed class (for the self-complementary C>G / T>A both do, so bp_is_ref is copied);
    otherwise copy peak SNP j's bp_is_ref (bit j = 1 means BP = ALT).
    """
    pool_ref = np.asarray(pool_ref)
    pool_alt = np.asarray(pool_alt)
    n_cand, kc = pool_ref.shape
    mask = np.zeros(n_cand, dtype=np.int64)
    for j in range(kc):
        pdc = DC_TABLE[peak_bp[j], peak_tp[j]]
        c0 = DC_TABLE[pool_ref[:, j], pool_alt[:, j]]          # pool class if BP = REF
        bit = np.full(n_cand, int(peak_bits[j]), dtype=np.int64)
        if not SELF_COMP[pdc]:
            match = UND_CLASS[c0] == UND_CLASS[pdc]
            bit = np.where(match, (c0 != pdc).astype(np.int64), bit)
        mask |= bit << j
    return mask


# --------------------------------------------------------------------------------------------------------------------
# Statistics helpers
# --------------------------------------------------------------------------------------------------------------------
def two_sided_p(obs, null):
    """Empirical two-sided p (README.md): min(1, 2 * min(lo, hi)) with lo = (#{null <= obs} + 1) / (N + 1) and
    hi = (#{null >= obs} + 1) / (N + 1). NaN null values (e.g. f with no events) are dropped from N.
    obs: scalar or (m,); null: (N,) or (N, m)."""
    obs = np.asarray(obs, dtype=float)
    null = np.asarray(null, dtype=float)
    scalar = obs.ndim == 0
    if null.ndim == 1:
        null = null[:, None]
    o = np.broadcast_to(np.atleast_1d(obs), (null.shape[1],))
    n = (~np.isnan(null)).sum(axis=0)
    lo = (null <= o).sum(axis=0)
    hi = (null >= o).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.minimum(1.0, 2.0 * np.minimum((lo + 1) / (n + 1), (hi + 1) / (n + 1)))
    p = np.where(np.isnan(o) | (n == 0), np.nan, p)
    return float(p[0]) if scalar else p


def bh_fdr(p):
    """Benjamini-Hochberg q-values (step-up, monotone, capped at 1); NaN entries are left out and stay NaN."""
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    m = int(ok.sum())
    if m == 0:
        return q
    pv = p[ok]
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


def frac(nbp, ntp):
    nbp = np.asarray(nbp, dtype=float)
    n = nbp + np.asarray(ntp, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, nbp / np.where(n > 0, n, 1), np.nan)


def ratio(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(b > 0, a / np.where(b > 0, b, 1), np.nan)


def describe(obs, null):
    null = np.asarray(null, dtype=float)
    v = null[~np.isnan(null)]
    if len(v) == 0:
        return dict(observed=float(obs), null_mean=np.nan, null_sd=np.nan, null_q025=np.nan, null_q975=np.nan,
                    n_valid_draws=0, p_two_sided=np.nan)
    return dict(observed=float(obs), null_mean=float(v.mean()), null_sd=float(v.std(ddof=1)) if len(v) > 1 else np.nan,
                null_q025=float(np.quantile(v, 0.025)), null_q975=float(np.quantile(v, 0.975)),
                n_valid_draws=int(len(v)), p_two_sided=two_sided_p(obs, null))


def sample_without_replacement(rng, n, g, n_draws):
    """(n_draws, g) indices into range(n): distinct within each row when g <= n (Floyd's algorithm, vectorised over
    draws, then a random column order), with replacement when g > n. Returns (idx, with_replacement)."""
    if n <= 0:
        raise ValueError("empty stratum")
    if g > n:
        return rng.integers(0, n, size=(n_draws, g)), True
    if g == 1:
        return rng.integers(0, n, size=(n_draws, 1)), False
    sel = np.empty((n_draws, g), dtype=np.int64)
    for c, j in enumerate(range(n - g, n)):
        t = rng.integers(0, j + 1, size=n_draws)
        if c:
            dup = (sel[:, :c] == t[:, None]).any(axis=1)
            t[dup] = j
        sel[:, c] = t
    perm = np.argsort(rng.random((n_draws, g)), axis=1)
    return np.take_along_axis(sel, perm, axis=1), False


def sd_guard(f_null, events_null, factor=SD_GUARD_FACTOR):
    """Null SD of f must be >= factor * 0.5 / sqrt(mean events) (binomial SD at f = 0.5)."""
    f = np.asarray(f_null, dtype=float)
    f = f[~np.isnan(f)]
    mean_n = float(np.mean(events_null)) if len(events_null) else 0.0
    if len(f) < 2 or mean_n <= 0:
        return dict(sd=np.nan, threshold=np.nan, mean_events=mean_n, passed=None)
    sd = float(np.std(f, ddof=1))
    thr = factor * 0.5 / np.sqrt(mean_n)
    return dict(sd=sd, threshold=float(thr), mean_events=mean_n, passed=bool(sd >= thr))


def calibration_pvalues(null, n_calib, rng):
    """Leave-one-out two-sided p of the first n_calib valid null draws against the other draws.
    Returns (p_raw (README.md formula), p_randomized (ties broken at random, exactly uniform under exchangeability))."""
    null = np.asarray(null, dtype=float)
    vals = null[~np.isnan(null)]
    if len(vals) < 3:
        return np.array([]), np.array([])
    srt = np.sort(vals)
    n_other = len(srt) - 1
    x = vals[:n_calib]
    le = np.searchsorted(srt, x, "right") - 1            # others <= x
    lt = np.searchsorted(srt, x, "left")                 # others < x
    ge = len(srt) - np.searchsorted(srt, x, "left") - 1  # others >= x
    ties = le - lt                                       # others == x
    p_raw = np.minimum(1.0, 2.0 * np.minimum((le + 1) / (n_other + 1), (ge + 1) / (n_other + 1)))
    p_lo = (lt + rng.random(len(x)) * (ties + 1)) / (n_other + 1)
    p_rand = 2.0 * np.minimum(p_lo, 1.0 - p_lo)
    return p_raw, p_rand


# --------------------------------------------------------------------------------------------------------------------
# Loading and validation
# --------------------------------------------------------------------------------------------------------------------
def _flat(series):
    return np.array(",".join(series.tolist()).split(","))


def load_units(path):
    """Read units.tsv and validate it against README.md. Returns a dict of per-unit and flattened per-SNP arrays."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_filter=False)
    missing = [c for c in UNIT_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing column(s) {missing}")
    if df.empty:
        raise ValueError(f"{path}: no units")
    try:
        uidx = df["unit_idx"].astype(np.int64)
    except ValueError as e:
        raise ValueError(f"{path}: unit_idx must be integers ({e})")
    dup = uidx[uidx.duplicated(keep=False)]
    if len(dup):
        raise ValueError(f"{path}: duplicate unit_idx ({dup.nunique()} value(s), e.g. {sorted(dup.unique())[:5]}); "
                         "every unit needs its own index (the collapsed-ID failure mode)")
    df["unit_idx"] = uidx
    df = df.sort_values("unit_idx").reset_index(drop=True)
    n = len(df)
    if not np.array_equal(df["unit_idx"].to_numpy(), np.arange(n)):
        raise ValueError(f"{path}: unit_idx must be 0..N-1 with no gaps")
    src = df["source"].to_numpy()
    if not np.isin(src, ["peak", "pool"]).all():
        raise ValueError(f"{path}: source must be peak or pool")
    is_peak = src == "peak"
    if is_peak.any() and (~is_peak).any() and np.where(is_peak)[0].max() > np.where(~is_peak)[0].min():
        warnings.warn("units.tsv: peak units are not all before pool units (specified order); continuing")
    if not ((df.loc[~is_peak, "locus"] == "pool").all() and (df.loc[~is_peak, "set"] == "pool").all()):
        raise ValueError(f"{path}: pool rows must have locus = set = pool")
    if not (df.loc[is_peak, "set"].isin(["extended", "core"]).all() and (df.loc[is_peak, "locus"] != "pool").all()):
        raise ValueError(f"{path}: peak rows need set extended/core and a locus name")
    k = df["k"].astype(np.int64).to_numpy()
    if (k < 1).any():
        raise ValueError(f"{path}: k must be >= 1")
    if (k > 15).any():
        raise ValueError(f"{path}: k up to 15 only (hits.npz masks are int16); max k = {k.max()}")
    for col in ["pos", "ref", "alt", "snp_off", "cpg", "af_ref", "snp_id"]:
        cnt = df[col].str.count(",").to_numpy() + 1
        if (cnt != k).any():
            i = int(np.where(cnt != k)[0][0])
            raise ValueError(f"{path}: column {col} has {cnt[i]} item(s) but k = {k[i]} (unit {i})")
    cnt = df.loc[is_peak, "log_p"].str.count(",").to_numpy() + 1
    if (cnt != k[is_peak]).any():
        raise ValueError(f"{path}: peak log_p must have k items")
    snp_ref = base_codes(_flat(df["ref"]))
    snp_alt = base_codes(_flat(df["alt"]))
    if (snp_ref == snp_alt).any():
        raise ValueError(f"{path}: REF equals ALT for some SNP")
    snp_pos = _flat(df["pos"]).astype(np.int64)
    snp_cpg = _flat(df["cpg"]).astype(np.int64)
    if not np.isin(snp_cpg, [0, 1]).all():
        raise ValueError(f"{path}: cpg must be 0/1")
    snp_af = _flat(df["af_ref"]).astype(float)
    snp_off0 = np.array([int(s.split(",", 1)[0]) for s in df["snp_off"]])
    if (snp_off0 != F + 1).any():
        raise ValueError(f"{path}: the first snp_off must be {F + 1}")
    snp_start = np.r_[0, np.cumsum(k)[:-1]].astype(np.int64)
    snp_unit = np.repeat(np.arange(n), k)
    # positions strictly increasing inside a unit, gaps <= F
    d = np.diff(snp_pos)
    same = snp_unit[1:] == snp_unit[:-1]
    if (d[same] <= 0).any() or (d[same] > F).any():
        raise ValueError(f"{path}: SNP positions in a unit must increase with gaps <= {F}")
    gc = pd.to_numeric(df["gc"], errors="coerce").to_numpy(dtype=float)
    if np.isnan(gc).any():
        raise ValueError(f"{path}: gc is missing for {int(np.isnan(gc).sum())} unit(s)")
    gap_str = df["gap_bin"].to_numpy()
    gap_bin = np.full(n, -1, dtype=np.int64)
    for i, g in enumerate(GAP_BINS):
        gap_bin[gap_str == g] = i
    if (gap_bin[k >= 2] < 0).any():
        raise ValueError(f"{path}: gap_bin must be one of {GAP_BINS} for k >= 2")
    if not np.isin(gap_str[k == 1], ["NA", ""]).all():
        raise ValueError(f"{path}: gap_bin must be NA for k = 1")
    bpm = pd.to_numeric(df["bp_mask"].replace({"NA": "", "nan": ""}), errors="coerce").to_numpy()
    if np.isnan(bpm[is_peak]).any():
        raise ValueError(f"{path}: peak rows need an integer bp_mask")
    if (~np.isnan(bpm[~is_peak])).any():
        raise ValueError(f"{path}: pool rows must have bp_mask = NA")
    bp_mask = np.where(is_peak, np.nan_to_num(bpm, nan=-1), -1).astype(np.int64)
    if ((bp_mask[is_peak] < 0) | (bp_mask[is_peak] >= (1 << k[is_peak]))).any():
        raise ValueError(f"{path}: bp_mask out of range 0..2^k-1")
    # a peak SNP may sit in only one unit per (locus, set); duplicated units would double-count calls
    pk = np.repeat(is_peak, k)
    key = pd.DataFrame({"locus": np.repeat(df["locus"].to_numpy(), k)[pk], "set": np.repeat(df["set"].to_numpy(), k)[pk],
                        "chr": np.repeat(df["chr"].to_numpy(), k)[pk], "pos": snp_pos[pk]})
    if key.duplicated().any():
        r = key[key.duplicated()].iloc[0]
        raise ValueError(f"{path}: peak SNP {r['chr']}:{r['pos']} appears in more than one {r['locus']}/{r['set']} unit")
    pool_key = (df.loc[~is_peak, "chr"] + ":" + df.loc[~is_peak, "pos"])
    if pool_key.duplicated().any():
        warnings.warn(f"units.tsv: {int(pool_key.duplicated().sum())} duplicated pool unit(s)")
    return dict(df=df, n=n, k=k, is_peak=is_peak, bp_mask=bp_mask, gc=gc, gap_bin=gap_bin, snp_start=snp_start,
                snp_unit=snp_unit, snp_ref=snp_ref, snp_alt=snp_alt, snp_pos=snp_pos, snp_cpg=snp_cpg, snp_af=snp_af,
                locus=df["locus"].to_numpy(), set=df["set"].to_numpy(), chr=df["chr"].to_numpy())


def load_hits(path, units):
    """Read hits.npz (README.md Stage F) and validate it against units."""
    try:
        z = np.load(path, allow_pickle=False)
        data = {key: z[key] for key in z.files}
    except ValueError:
        z = np.load(path, allow_pickle=True)   # string arrays saved with dtype=object
        data = {key: z[key] for key in z.files}
    for key in ["unit_idx", "mask", "motif_idx", "mlog10p", "motif_ids"]:
        if key not in data:
            raise ValueError(f"{path}: missing array {key}")
    u = data["unit_idx"].astype(np.int64)
    m = data["mask"].astype(np.int64)
    mot = data["motif_idx"].astype(np.int64)
    s = data["mlog10p"].astype(np.float64)
    if not (len(u) == len(m) == len(mot) == len(s)):
        raise ValueError(f"{path}: unit_idx, mask, motif_idx and mlog10p differ in length")
    motif_ids = np.array([str(x) for x in data["motif_ids"]])
    motif_names = np.array([str(x) for x in data["motif_names"]]) if "motif_names" in data else motif_ids.copy()
    if len(set(motif_ids.tolist())) != len(motif_ids):
        raise ValueError(f"{path}: motif_ids are not unique")
    if len(u):
        if u.min() < 0 or u.max() >= units["n"]:
            raise ValueError(f"{path}: unit_idx outside 0..{units['n'] - 1}; hits and units.tsv are out of sync")
        if mot.min() < 0 or mot.max() >= len(motif_ids):
            raise ValueError(f"{path}: motif_idx outside 0..{len(motif_ids) - 1}")
        kk = units["k"][u]
        if ((m < 0) | (m >= (1 << kk))).any():
            raise ValueError(f"{path}: mask outside 0..2^k-1 for some hit")
        if not np.isfinite(s).all():
            raise ValueError(f"{path}: non-finite mlog10p")
    low = s < SCAN_MLOG10P - MLOG10P_TOL
    if low.any():
        warnings.warn(f"hits.npz: {int(low.sum())} entries with p > 1e-3 ignored (scan threshold is 1e-3)")
        keep = ~low
        u, m, mot, s = u[keep], m[keep], mot[keep], s[keep]
    return dict(unit_idx=u, mask=m, motif_idx=mot, mlog10p=s, motif_ids=motif_ids, motif_names=motif_names)


def _primary(x):
    x = str(x).split(";")[0].strip()
    return x if x and x.lower() not in ("nan", "na", "none") else "Unknown"


def load_families(path, motif_ids, motif_names):
    fam = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_filter=False)
    for c in ["motif_id", "motif_name", "tf_class", "tf_family"]:
        if c not in fam.columns:
            raise ValueError(f"{path}: missing column {c}")
    if fam["motif_id"].duplicated().any():
        raise ValueError(f"{path}: duplicated motif_id {fam.loc[fam['motif_id'].duplicated(), 'motif_id'].iloc[0]}")
    fam = fam.set_index("motif_id")
    missing = [x for x in motif_ids if x not in fam.index]
    if missing:
        raise ValueError(f"{path}: {len(missing)} motif(s) of hits.npz have no family row, e.g. {missing[:3]}")
    sub = fam.loc[list(motif_ids)]
    primary = np.array([_primary(x) for x in sub["tf_family"]])
    pclass = np.array([_primary(x) for x in sub["tf_class"]])
    all_fam = np.array([x if x.strip() else "Unknown" for x in sub["tf_family"]])
    names = sorted(set(primary.tolist()))
    motif_fam = np.array([names.index(x) for x in primary], dtype=np.int64)
    fam_class = []
    for i, nm in enumerate(names):
        cl = pd.Series(pclass[motif_fam == i]).value_counts()
        fam_class.append(cl.index[0] if len(cl) else "Unknown")
    return dict(names=names, n=len(names), motif_fam=motif_fam, primary=primary, all_families=all_fam,
                motif_class=pclass, fam_class=fam_class, fam_n_motifs=np.bincount(motif_fam, minlength=len(names)))


def named_group_map(motif_ids):
    names = list(NAMED_GROUPS)
    motif_grp = np.full(len(motif_ids), -1, dtype=np.int64)
    found = {}
    index = {x: i for i, x in enumerate(motif_ids)}
    for g, (name, members) in enumerate(NAMED_GROUPS.items()):
        found[name] = [x for x in members if x in index]
        lost = [x for x in members if x not in index]
        if lost:
            warnings.warn(f"named group {name}: motif(s) {lost} not in hits.npz motif_ids; tested without them")
        for x in found[name]:
            if motif_grp[index[x]] >= 0:
                raise ValueError(f"motif {x} is in two named groups")
            motif_grp[index[x]] = g
    return dict(names=names, n=len(names), motif_grp=motif_grp, found=found)


def load_polarized(path):
    pol = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_filter=False)
    for c in ["chr", "pos", "ref", "alt", "bp_is_ref"]:
        if c not in pol.columns:
            raise ValueError(f"{path}: missing column {c}")
    pol["pos"] = pol["pos"].astype(np.int64)
    pol["bp_is_ref_b"] = pol["bp_is_ref"].str.strip().str.lower().isin(["true", "1", "t", "yes"])
    return pol


def check_polarized(units, pol):
    """bp_mask must agree with peak_snps_polarized.tsv bp_is_ref for every peak SNP (orientation is never inferred
    from REF)."""
    k = units["k"]
    pk = np.repeat(units["is_peak"], k)
    j = np.arange(len(units["snp_unit"])) - units["snp_start"][units["snp_unit"]]
    bit = (units["bp_mask"][units["snp_unit"]] >> j) & 1
    peak = pd.DataFrame({"chr": units["chr"][units["snp_unit"]][pk], "pos": units["snp_pos"][pk],
                         "ref": np.array(list(BASES))[units["snp_ref"][pk]],
                         "alt": np.array(list(BASES))[units["snp_alt"][pk]], "bp_is_ref_u": bit[pk] == 0})
    mrg = peak.merge(pol.drop_duplicates(["chr", "pos"]), on=["chr", "pos"], how="left", suffixes=("", "_pol"),
                     indicator=True)
    n_miss = int((mrg["_merge"] == "left_only").sum())
    if n_miss:
        warnings.warn(f"{n_miss} peak SNP(s) of units.tsv are not in peak_snps_polarized.tsv")
    both = mrg[mrg["_merge"] == "both"]
    bad_allele = both[(both["ref"] != both["ref_pol"].str.upper()) | (both["alt"] != both["alt_pol"].str.upper())]
    if len(bad_allele):
        r = bad_allele.iloc[0]
        raise ValueError(f"REF/ALT differ between units.tsv and peak_snps_polarized.tsv at {r['chr']}:{r['pos']}")
    bad = both[both["bp_is_ref_u"] != both["bp_is_ref_b"]]
    if len(bad):
        r = bad.iloc[0]
        raise ValueError(f"bp_mask disagrees with peak_snps_polarized.tsv bp_is_ref at {len(bad)} SNP(s), "
                         f"e.g. {r['chr']}:{r['pos']}")
    return len(both), n_miss


# --------------------------------------------------------------------------------------------------------------------
# Haplotype rows and directional tables
# --------------------------------------------------------------------------------------------------------------------
def build_rows(units):
    """One row per haplotype: pool units get 2^k rows (row = mask), peak units 2 rows (0 = bp_mask, 1 = complement).
    A row read as 'BP haplotype' has its TP haplotype at flip[row]."""
    k = units["k"]
    is_peak = units["is_peak"]
    nslots = np.where(is_peak, 2, 1 << k).astype(np.int64)
    hap_off = np.r_[0, np.cumsum(nslots)[:-1]].astype(np.int64)
    n = int(nslots.sum())
    row_unit = np.repeat(np.arange(units["n"]), nslots)
    row_slot = np.arange(n) - hap_off[row_unit]
    full = (1 << k[row_unit]) - 1
    rp = is_peak[row_unit]
    bpm = units["bp_mask"][row_unit]
    row_mask = np.where(rp, np.where(row_slot == 0, bpm, bpm ^ full), row_slot)
    flip = np.where(rp, hap_off[row_unit] + 1 - row_slot, hap_off[row_unit] + (row_slot ^ full))
    return dict(n=n, hap_off=hap_off, row_unit=row_unit, row_mask=row_mask, flip=flip)


def hits_to_rows(hits, units, rows):
    u, m = hits["unit_idx"], hits["mask"]
    peak = units["is_peak"][u]
    bpm = units["bp_mask"][u]
    full = (1 << units["k"][u]) - 1
    bad = peak & (m != bpm) & (m != (bpm ^ full))
    if bad.any():
        i = int(np.where(bad)[0][0])
        raise ValueError(f"hits.npz: {int(bad.sum())} peak hit(s) on masks other than bp_mask or its complement "
                         f"(e.g. unit {u[i]} mask {m[i]}); hits.npz and units.tsv are out of sync")
    slot = np.where(peak, np.where(m == bpm, 0, 1), m)
    hap = rows["hap_off"][u] + slot
    n_mot = len(hits["motif_ids"])
    key = np.sort(hap * n_mot + hits["motif_idx"])
    if len(key) > 1 and (np.diff(key) == 0).any():
        i = int(np.where(np.diff(key) == 0)[0][0])
        raise ValueError(f"hits.npz: duplicated (unit, mask, motif) entry (row key {key[i]}); the format stores one "
                         "minimum-p entry per triple")
    return hap


def _isin_sorted(q, srt):
    if len(srt) == 0 or len(q) == 0:
        return np.zeros(len(q), dtype=bool)
    i = np.searchsorted(srt, q)
    i = np.minimum(i, len(srt) - 1)
    return srt[i] == q


def _category_tables(p_hap, p_cat, w_hap, w_cat, n_cat, flip, n_rows):
    """Collapse motif hits to categories (family / named group). Returns per rule: B (rows x cats, 1 = category is
    BP-only when the row is read as BP) and Bt = B[flip] (TP-only)."""
    out = {}
    pk = np.unique(p_hap * n_cat + p_cat)
    wk = np.unique(w_hap * n_cat + w_cat)
    h, c = pk // n_cat, pk % n_cat
    fk = flip[h] * n_cat + c
    for rule, other in (("standard", pk), ("strong", wk)):
        keep = ~_isin_sorted(fk, other)
        B = sparse.csr_matrix((np.ones(int(keep.sum()), dtype=np.int32), (h[keep], c[keep])), shape=(n_rows, n_cat))
        out[rule] = dict(B=B, Bt=B[flip].tocsr())
    return out


def directional_tables(hap, motif, mlog10p, rows, motif_fam, n_fam, motif_grp, n_grp, n_mot):
    """Per haplotype row r read as BP (TP = flip[r]): motif-level BP-only counts for both rules and the family /
    named-group category matrices."""
    flip = rows["flip"]
    n_rows = rows["n"]
    anyhit = mlog10p >= SCAN_MLOG10P - MLOG10P_TOL
    pres = mlog10p > PRESENT_MLOG10P
    w_hap, w_mot = hap[anyhit], motif[anyhit]
    p_hap, p_mot = hap[pres], motif[pres]
    p_keys = np.sort(p_hap * n_mot + p_mot)
    w_keys = np.sort(w_hap * n_mot + w_mot)
    fk = flip[p_hap] * n_mot + p_mot
    in_pt = _isin_sorted(fk, p_keys)
    in_wt = _isin_sorted(fk, w_keys)
    T = {}
    for rule, other in (("standard", in_pt), ("strong", in_wt)):
        nbp = np.bincount(p_hap[~other], minlength=n_rows).astype(np.int64)
        T[rule] = dict(nBP=nbp, nTP=nbp[flip])
    fam = _category_tables(p_hap, motif_fam[p_mot], w_hap, motif_fam[w_mot], n_fam, flip, n_rows)
    gp, gw = motif_grp[p_mot] >= 0, motif_grp[w_mot] >= 0
    grp = _category_tables(p_hap[gp], motif_grp[p_mot[gp]], w_hap[gw], motif_grp[w_mot[gw]], n_grp, flip, n_rows)
    for rule in RULES:
        T[rule]["fam"] = fam[rule]
        T[rule]["grp"] = grp[rule]
        T[rule]["nBPfam"] = np.diff(fam[rule]["B"].indptr).astype(np.int64)
        T[rule]["nTPfam"] = T[rule]["nBPfam"][flip]
    return T


def row_covariates(units, rows):
    """Balance covariates per row read as BP: k, gc and per-SNP sums of cpg, af_ref, bp_is_ref and class counts."""
    ru = rows["row_unit"]
    kr = units["k"][ru]
    tot = int(kr.sum())
    rep = np.repeat(np.arange(rows["n"]), kr)
    j = np.arange(tot) - np.repeat(np.r_[0, np.cumsum(kr)[:-1]], kr)
    snp = units["snp_start"][ru[rep]] + j
    bit = (rows["row_mask"][rep] >> j) & 1
    ref, alt = units["snp_ref"][snp], units["snp_alt"][snp]
    cls = DC_TABLE[np.where(bit == 1, alt, ref), np.where(bit == 1, ref, alt)]
    n = rows["n"]
    return dict(k=kr.astype(np.int64), gc=units["gc"][ru],
                cpg=np.bincount(rep, weights=units["snp_cpg"][snp], minlength=n),
                af=np.bincount(rep, weights=units["snp_af"][snp], minlength=n),
                bpref=np.bincount(rep, weights=(bit == 0).astype(float), minlength=n),
                cls=np.bincount(rep * 6 + cls, minlength=n * 6).reshape(n, 6).astype(float))


def peak_bases(units, u):
    """BP base, TP base and bp_mask bit per SNP of peak unit u (orientation from bp_mask only)."""
    s0, k = units["snp_start"][u], units["k"][u]
    ref = units["snp_ref"][s0:s0 + k]
    alt = units["snp_alt"][s0:s0 + k]
    bits = (units["bp_mask"][u] >> np.arange(k)) & 1
    return np.where(bits == 1, alt, ref), np.where(bits == 1, ref, alt), bits


# --------------------------------------------------------------------------------------------------------------------
# Strata
# --------------------------------------------------------------------------------------------------------------------
class StratumIndex:
    """Members (ascending ids) indexed by integer key fields at several fallback levels."""

    def __init__(self, member_ids, fields, radix, levels):
        self.member_ids = np.asarray(member_ids)
        self.radix = radix
        self.levels = levels
        self.tables = []
        for _, keys in levels:
            code = self._encode(fields, keys)
            order = np.argsort(code, kind="stable")
            uniq, start, cnt = np.unique(code[order], return_index=True, return_counts=True)
            self.tables.append((keys, uniq, start, cnt, self.member_ids[order]))

    def _encode(self, vals, keys):
        code = 0
        for key in keys:
            code = code * self.radix[key] + np.asarray(vals[key], dtype=np.int64)
        return code

    def members(self, level, query):
        keys, uniq, start, cnt, mem = self.tables[level]
        c = int(self._encode(query, keys))
        i = int(np.searchsorted(uniq, c))
        if i < len(uniq) and uniq[i] == c:
            return mem[start[i]:start[i] + cnt[i]], c
        return mem[:0], c

    def lookup(self, query, min_size):
        """First level with >= min_size members (af -> gc -> cpg order is the order of self.levels). If none,
        the coarsest non-empty level, flagged below_min. Returns (level, members, code, sizes, below_min)."""
        sizes = []
        for lev in range(len(self.levels)):
            mem, code = self.members(lev, query)
            sizes.append(len(mem))
            if len(mem) >= min_size:
                return lev, mem, code, sizes, False
        if len(mem) == 0:
            raise ValueError(f"no pool units at all for stratum {query}")
        return len(self.levels) - 1, mem, code, sizes, True


K1_LEVELS = [("full", ["cls", "cpg", "gcq", "afq", "bpref"]),
             ("drop_af", ["cls", "cpg", "gcq", "bpref"]),
             ("drop_af_gc", ["cls", "cpg", "bpref"]),
             ("drop_af_gc_cpg", ["cls", "bpref"])]
K2_LEVELS = [("full", ["kc", "gap", "gcq"]),
             ("drop_gc", ["kc", "gap"]),
             ("drop_gc_gap", ["kc"])]
RADIX = dict(cls=6, cpg=2, gcq=5, afq=5, bpref=2, kc=16, gap=3)


def quintile_edges(x):
    return np.quantile(np.asarray(x, dtype=float), [0.2, 0.4, 0.6, 0.8])


def qbin(x, edges):
    return np.searchsorted(edges, x, side="right")


def key_string(level_keys, query):
    parts = []
    for key in level_keys:
        v = int(query[key])
        parts.append(f"class={CLASSES[v]}" if key == "cls" else f"gap_bin={GAP_BINS[v]}" if key == "gap"
                     else f"bp_is_ref={v}" if key == "bpref" else f"{key}={v}")
    return ";".join(parts)


def build_strata(units, rows):
    pool = ~units["is_peak"]
    p1 = np.where(pool & (units["k"] == 1))[0]
    p2 = np.where(pool & (units["k"] >= 2))[0]
    if len(p1) == 0:
        raise ValueError("units.tsv has no pool singletons")
    snp1 = units["snp_start"][p1]
    gc1_edges = quintile_edges(units["gc"][p1])
    af_edges = quintile_edges(units["snp_af"][snp1])
    # every pool singleton in both orientations: row hap_off+0 is BP = REF, hap_off+1 is BP = ALT
    r = np.sort(np.r_[rows["hap_off"][p1], rows["hap_off"][p1] + 1])
    u = rows["row_unit"][r]
    bit = rows["row_mask"][r]
    s = units["snp_start"][u]
    ref, alt = units["snp_ref"][s], units["snp_alt"][s]
    f1 = dict(cls=DC_TABLE[np.where(bit == 1, alt, ref), np.where(bit == 1, ref, alt)], cpg=units["snp_cpg"][s],
              gcq=qbin(units["gc"][u], gc1_edges), afq=qbin(units["snp_af"][s], af_edges), bpref=1 - bit)
    idx1 = StratumIndex(r, f1, RADIX, K1_LEVELS)
    out = dict(idx1=idx1, gc1_edges=gc1_edges, af_edges=af_edges, idx2=None, gc2_edges=None, kmax=1,
               n_pool1=len(p1), n_pool2=len(p2), pool2_k=np.bincount(units["k"][p2]) if len(p2) else np.array([]))
    if len(p2):
        gc2_edges = quintile_edges(units["gc"][p2])
        f2 = dict(kc=units["k"][p2], gap=units["gap_bin"][p2], gcq=qbin(units["gc"][p2], gc2_edges))
        out.update(idx2=StratumIndex(p2, f2, RADIX, K2_LEVELS), gc2_edges=gc2_edges, kmax=int(units["k"][p2].max()))
    return out


def make_draws(units, rows, strata, peak_u, n_draws, rng, min_stratum, locus, set_name):
    """Matched null draws for one locus and set: (n_draws, n_peak_units) BP-row indices + one log row per peak unit.

    Groups of peak units that share a final stratum are sampled together without replacement inside a draw. Groups
    are processed in order of their first peak unit and members are in ascending id order, so swapping BP and TP at
    the peak mirrors every draw exactly (same pool units, flipped orientation)."""
    plans, log_rows = [], []
    kmax = strata["kmax"]
    for u in peak_u:
        u = int(u)
        k = int(units["k"][u])
        bp, tp, bits = peak_bases(units, u)
        s0 = units["snp_start"][u]
        if k == 1:
            q = dict(cls=int(DC_TABLE[bp[0], tp[0]]), cpg=int(units["snp_cpg"][s0]),
                     gcq=int(qbin(units["gc"][u], strata["gc1_edges"])),
                     afq=int(qbin(units["snp_af"][s0], strata["af_edges"])), bpref=int(1 - bits[0]))
            lev, mem, code, sizes, below = strata["idx1"].lookup(q, min_stratum)
            cand = mem
            levels = K1_LEVELS
            gkey = ("k1", lev, code)
        else:
            if strata["idx2"] is None:
                raise ValueError(f"peak unit {u} has k = {k} but units.tsv has no pool chains")
            kc = min(k, kmax)
            q = dict(kc=kc, gap=int(units["gap_bin"][u]), gcq=int(qbin(units["gc"][u], strata["gc2_edges"])))
            lev, mem, code, sizes, below = strata["idx2"].lookup(q, min_stratum)
            kc_pool = int(units["k"][mem[0]])
            snp = units["snp_start"][mem][:, None] + np.arange(kc_pool)[None, :]
            masks = orient_chain_masks(units["snp_ref"][snp], units["snp_alt"][snp], bp, tp, bits)
            cand = rows["hap_off"][mem] + masks
            levels = K2_LEVELS
            gkey = ("k2", lev, code)
        plans.append(cand)
        log_rows.append(dict(locus=locus, set=set_name, unit_idx=u, snp_id=units["df"].at[u, "snp_id"], k=k,
                             pool_k=1 if k == 1 else int(units["k"][mem[0]]), k_capped=k > 1 and k > kmax,
                             stratum_type="k1" if k == 1 else "k2", level=levels[lev][0], fallback=lev > 0,
                             below_min=below, stratum_key=key_string(levels[lev][1], q), pool_size=len(cand),
                             sizes_by_level="/".join(str(x) for x in sizes), _gkey=gkey))
    or_idx = np.empty((n_draws, len(peak_u)), dtype=np.int64)
    groups = {}
    for i, row in enumerate(log_rows):
        groups.setdefault(row["_gkey"], []).append(i)
    for members in groups.values():
        n = len(plans[members[0]])
        pos, repl = sample_without_replacement(rng, n, len(members), n_draws)
        for slot, i in enumerate(members):
            or_idx[:, i] = plans[i][pos[:, slot]]
            log_rows[i]["group_size"] = len(members)
            log_rows[i]["with_replacement"] = repl
    for row in log_rows:
        del row["_gkey"]
    return or_idx, log_rows


# --------------------------------------------------------------------------------------------------------------------
# Aggregation
# --------------------------------------------------------------------------------------------------------------------
def aggregate(or_idx, T, cov, n_rows, n_fam, n_grp):
    """Per-draw totals for draws of BP rows (n_draws, n_units)."""
    nd, nu = or_idx.shape
    D = sparse.csr_matrix((np.ones(nd * nu, dtype=np.int32), (np.repeat(np.arange(nd), nu), or_idx.ravel())),
                          shape=(nd, n_rows))
    a = dict(n_units=nu, nsnp=cov["k"][or_idx].sum(1), sum_k=cov["k"][or_idx].sum(1).astype(float),
             sum_gc=cov["gc"][or_idx].sum(1), sum_cpg=cov["cpg"][or_idx].sum(1), sum_af=cov["af"][or_idx].sum(1),
             sum_bpref=cov["bpref"][or_idx].sum(1), sum_cls=np.asarray(D @ cov["cls"]))
    for rule in RULES:
        t = T[rule]
        a[rule] = dict(nBP=t["nBP"][or_idx].sum(1), nTP=t["nTP"][or_idx].sum(1),
                       nBPfam=t["nBPfam"][or_idx].sum(1), nTPfam=t["nTPfam"][or_idx].sum(1),
                       famBP=(D @ t["fam"]["B"]).toarray(), famTP=(D @ t["fam"]["Bt"]).toarray(),
                       grpBP=(D @ t["grp"]["B"]).toarray(), grpTP=(D @ t["grp"]["Bt"]).toarray())
    return a


def add_aggregates(aggs):
    def _add(xs):
        if isinstance(xs[0], dict):
            return {key: _add([x[key] for x in xs]) for key in xs[0]}
        return sum(xs[1:], xs[0])
    return _add(aggs)


# --------------------------------------------------------------------------------------------------------------------
# Main analysis
# --------------------------------------------------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(description="TFBS Stage B: matched-null test (README.md).")
    p.add_argument("--units", required=True, help="$OUT/units.tsv")
    p.add_argument("--hits", required=True, help="$OUT/hits.npz")
    p.add_argument("--families", required=True, help="input_data/15_TFBS_Scan/jaspar2024_families.tsv")
    p.add_argument("--polarized", default=None, help="$OUT/peak_snps_polarized.tsv (optional cross-check)")
    p.add_argument("--out_dir", required=True, help="$OUT/results")
    p.add_argument("--n_draws", type=int, default=10000)
    p.add_argument("--n_calib", type=int, default=500)
    p.add_argument("--seed", type=int, default=20261013)
    p.add_argument("--min_stratum", type=int, default=50)
    p.add_argument("--n_power", type=int, default=100, help="pseudo-peaks for the power estimate")
    p.add_argument("--no_plots", action="store_true")
    a = p.parse_args(argv)
    if a.n_calib >= a.n_draws:
        p.error("--n_calib must be smaller than --n_draws")
    return a


def _rng(seed, *key):
    return np.random.default_rng(np.random.SeedSequence(seed, spawn_key=tuple(int(x) for x in key)))


def run(args):
    t0 = time.time()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    log(f"loading {args.units}")
    units = load_units(args.units)
    log(f"loading {args.hits}")
    hits = load_hits(args.hits, units)
    fam = load_families(args.families, hits["motif_ids"], hits["motif_names"])
    grp = named_group_map(hits["motif_ids"])
    pol, pol_check = None, None
    if args.polarized:
        pol = load_polarized(args.polarized)
        pol_check = check_polarized(units, pol)
        log(f"bp_mask agrees with {args.polarized} at {pol_check[0]} SNPs ({pol_check[1]} not found)")
    rows = build_rows(units)
    hap = hits_to_rows(hits, units, rows)
    n_mot = len(hits["motif_ids"])
    log(f"{units['n']} units, {rows['n']} haplotype rows, {len(hap)} hits, {n_mot} motifs, {fam['n']} families")
    T = directional_tables(hap, hits["motif_idx"], hits["mlog10p"], rows, fam["motif_fam"], fam["n"],
                           grp["motif_grp"], grp["n"], n_mot)
    cov = row_covariates(units, rows)
    strata = build_strata(units, rows)
    log(f"directional tables + strata built ({time.time() - t0:.1f}s); KMAX = {strata['kmax']}")

    peak = units["is_peak"]
    present = set(units["locus"][peak].tolist())
    loci = [x for x in LOCUS_ORDER if x in present] + sorted(present - set(LOCUS_ORDER))
    sets = [s for s in SETS if (units["set"][peak] == s).any()]
    if "extended" not in sets:
        raise ValueError("units.tsv has no extended peak units")
    max_peak_k = int(units["k"][peak].max())
    if max_peak_k >= 2 and strata["kmax"] != min(4, max_peak_k):
        warnings.warn(f"pool KMAX = {strata['kmax']} but specified KMAX = min(4, largest peak k) = {min(4, max_peak_k)}")
    peak_u = {(s, L): np.where(peak & (units["set"] == s) & (units["locus"] == L))[0] for s in sets for L in loci}

    # ---- draws and aggregation
    draws, fb_rows = {}, []
    for si, s in enumerate(sets):
        for li, L in enumerate(loci):
            rng = _rng(args.seed, 1, si, li)
            draws[(s, L)], lr = make_draws(units, rows, strata, peak_u[(s, L)], args.n_draws, rng, args.min_stratum, L, s)
            fb_rows += lr
    fallbacks = pd.DataFrame(fb_rows)
    fallbacks.to_csv(out / "fallbacks.tsv", sep="\t", index=False)
    log(f"draws done ({time.time() - t0:.1f}s)")
    A = {}
    for s in sets:
        for L in loci:
            A[(s, L, "peak")] = aggregate(rows["hap_off"][peak_u[(s, L)]][None, :], T, cov, rows["n"], fam["n"], grp["n"])
            A[(s, L, "null")] = aggregate(draws[(s, L)], T, cov, rows["n"], fam["n"], grp["n"])
        for kind in ("peak", "null"):
            A[(s, "pooled", kind)] = add_aggregates([A[(s, L, kind)] for L in loci])
    scopes = loci + ["pooled"]
    log(f"aggregation done ({time.time() - t0:.1f}s)")

    # ---- sign-flip null (per set and rule; one flip matrix per set shared by loci and pooled)
    flips = {}
    for si, s in enumerate(sets):
        allu = np.concatenate([peak_u[(s, L)] for L in loci])
        lab = np.concatenate([[L] * len(peak_u[(s, L)]) for L in loci])
        r = rows["hap_off"][allu]
        flip = _rng(args.seed, 2, si).random((args.n_draws, len(allu))) < 0.5
        sign = 1.0 - 2.0 * flip
        for rule in RULES:
            t = T[rule]
            dlt = (t["nTP"][r] - t["nBP"][r]).astype(float)
            fam_s = (t["fam"]["B"][r] - t["fam"]["Bt"][r]).toarray().astype(float)
            grp_s = (t["grp"]["B"][r] - t["grp"]["Bt"][r]).toarray().astype(float)
            for sc in scopes:
                cols = np.ones(len(allu), bool) if sc == "pooled" else lab == sc
                nbp = t["nBP"][r][cols].sum()
                ntot = nbp + t["nTP"][r][cols].sum()
                fl = flip[:, cols].astype(float)
                flips[(s, rule, sc)] = dict(
                    f=frac(nbp + fl @ dlt[cols], ntot - nbp - fl @ dlt[cols]),
                    famnet=sign[:, cols] @ fam_s[cols], grpnet=sign[:, cols] @ grp_s[cols])

    # ---- locus tests
    lt_rows, null_f = [], {}
    for s in sets:
        for rule in RULES:
            for sc in scopes:
                pk, nl = A[(s, sc, "peak")], A[(s, sc, "null")]
                pr, nr = pk[rule], nl[rule]
                f_obs = frac(pr["nBP"], pr["nTP"])[0]
                f_null = frac(nr["nBP"], nr["nTP"])
                null_f[(s, rule, sc)] = (f_null, nr["nBP"] + nr["nTP"])
                base = dict(set=s, rule=rule, scope=sc, peak_units=pk["n_units"], peak_snps=int(pk["nsnp"][0]),
                            peak_nBP=int(pr["nBP"][0]), peak_nTP=int(pr["nTP"][0]),
                            null_mean_nBP=float(nr["nBP"].mean()), null_mean_nTP=float(nr["nTP"].mean()))
                lt_rows.append(dict(base, stat="f", null="matched", **describe(f_obs, f_null)))
                lt_rows.append(dict(base, stat="events_per_snp", null="matched",
                                    **describe(ratio(pr["nBP"] + pr["nTP"], pk["nsnp"])[0],
                                               ratio(nr["nBP"] + nr["nTP"], nl["nsnp"]))))
                lt_rows.append(dict(base, stat="f_family", null="matched",
                                    **describe(frac(pr["nBPfam"], pr["nTPfam"])[0], frac(nr["nBPfam"], nr["nTPfam"]))))
                lt_rows.append(dict(base, stat="f", null="sign_flip", **describe(f_obs, flips[(s, rule, sc)]["f"])))
            # Equal-weight pooled f: mean of the per-locus f, so the 200-kb cntnap5 window does not dominate.
            # Draw d at each locus forms joint draw d (README.md), so the null averages aligned draws.
            if len(loci) >= 2:
                fe_obs = float(np.mean([frac(A[(s, L, "peak")][rule]["nBP"], A[(s, L, "peak")][rule]["nTP"])[0]
                                        for L in loci]))
                fe_null = np.nanmean(np.vstack([null_f[(s, rule, L)][0] for L in loci]), axis=0)
                pk = A[(s, "pooled", "peak")]
                lt_rows.append(dict(set=s, rule=rule, scope="pooled", stat="f_equal_weight", null="matched",
                                    peak_units=pk["n_units"], peak_snps=int(pk["nsnp"][0]),
                                    peak_nBP=int(pk[rule]["nBP"][0]), peak_nTP=int(pk[rule]["nTP"][0]),
                                    **describe(fe_obs, fe_null)))
            if len(loci) >= 2:
                keepf = np.array([x != "Unknown" for x in fam["names"]])
                obs_pos = sum((A[(s, L, "peak")][rule]["famBP"] - A[(s, L, "peak")][rule]["famTP"] > 0).astype(int)
                              for L in loci)
                nul_pos = sum((A[(s, L, "null")][rule]["famBP"] - A[(s, L, "null")][rule]["famTP"] > 0).astype(int)
                              for L in loci)
                r_obs = int((obs_pos[0, keepf] >= 2).sum())
                r_null = (nul_pos[:, keepf] >= 2).sum(1)
                rec_fams = [fam["names"][i] for i in np.where(keepf & (obs_pos[0] >= 2))[0]]
                lt_rows.append(dict(set=s, rule=rule, scope="all_loci", stat="recurrence", null="matched_joint",
                                    detail=";".join(rec_fams), **describe(r_obs, r_null)))
    locus_tests = pd.DataFrame(lt_rows)
    lead = ["set", "rule", "scope", "stat", "null", "observed", "null_mean", "null_sd", "null_q025", "null_q975",
            "n_valid_draws", "p_two_sided"]
    locus_tests = locus_tests[lead + [c for c in locus_tests.columns if c not in lead]]
    locus_tests.to_csv(out / "locus_tests.tsv", sep="\t", index=False)

    # ---- family tests
    ft_rows, fam_null = [], {}
    for s in sets:
        for rule in RULES:
            for sc in scopes:
                pk, nl = A[(s, sc, "peak")][rule], A[(s, sc, "null")][rule]
                obs = (pk["famBP"] - pk["famTP"])[0]
                nul = nl["famBP"] - nl["famTP"]
                fam_null[(s, rule, sc)] = (obs, nul)
                testable = (obs != 0) | (nul != 0).any(0)
                p = np.where(testable, two_sided_p(obs, nul), np.nan)
                fl = flips[(s, rule, sc)]["famnet"]
                p_fl = np.where(testable & (obs != 0), two_sided_p(obs, fl), np.nan)
                q, q_fl = bh_fdr(p), bh_fdr(p_fl)
                sd = nul.std(0, ddof=1)
                for i, nm in enumerate(fam["names"]):
                    if not testable[i]:
                        continue
                    ft_rows.append(dict(set=s, rule=rule, scope=sc, tf_family=nm, tf_class=fam["fam_class"][i],
                                        n_motifs=int(fam["fam_n_motifs"][i]), peak_nBP=int(pk["famBP"][0, i]),
                                        peak_nTP=int(pk["famTP"][0, i]), peak_net=int(obs[i]),
                                        null_mean_net=float(nul[:, i].mean()), null_sd_net=float(sd[i]),
                                        null_q025=float(np.quantile(nul[:, i], 0.025)),
                                        null_q975=float(np.quantile(nul[:, i], 0.975)),
                                        p_two_sided=p[i], q_bh=q[i], signflip_p=p_fl[i], signflip_q=q_fl[i]))
    family_tests = pd.DataFrame(ft_rows)
    if len(family_tests):
        family_tests = family_tests.sort_values(["set", "rule", "scope", "q_bh", "p_two_sided", "tf_family"])
    family_tests.to_csv(out / "family_tests.tsv", sep="\t", index=False)

    # ---- per-unit calls (peak, motif level) + named-group member counts
    calls = peak_calls(units, hits, hap, rows, fam, grp, pol)
    for rule in RULES:   # internal consistency: the call table reproduces the per-row counts
        col = "direction" if rule == "standard" else "strong_direction"
        for s in sets:
            for L in loci:
                sub = calls[(calls["set"] == s) & (calls["locus"] == L)]
                pr = A[(s, L, "peak")][rule]
                if (int((sub[col] == "BP-only").sum()) != int(pr["nBP"][0])
                        or int((sub[col] == "TP-only").sum()) != int(pr["nTP"][0])):
                    raise RuntimeError(f"internal check failed: per-unit calls disagree with counts ({s}/{L}/{rule})")
    for L in loci:
        sub = calls[calls["locus"] == L].drop(columns=["_unit_order"])
        sub.to_csv(out / f"per_unit_calls_{L}.tsv", sep="\t", index=False)

    # ---- named groups (statistic 4b): Bonferroni across groups inside each scope; outside the family BH
    ng_rows = []
    for s in sets:
        for rule in RULES:
            for sc in scopes:
                pk, nl = A[(s, sc, "peak")], A[(s, sc, "null")]
                pr, nr = pk[rule], nl[rule]
                fl = flips[(s, rule, sc)]["grpnet"]
                for g, gname in enumerate(grp["names"]):
                    members = ",".join(grp["found"][gname])
                    net_o = int(pr["grpBP"][0, g] - pr["grpTP"][0, g])
                    net_n = nr["grpBP"][:, g] - nr["grpTP"][:, g]
                    eps_o = ratio(pr["grpBP"][0, g] + pr["grpTP"][0, g], pk["nsnp"][0])
                    eps_n = ratio(nr["grpBP"][:, g] + nr["grpTP"][:, g], nl["nsnp"])
                    base = dict(set=s, rule=rule, scope=sc, group=gname, members=members,
                                peak_nBP=int(pr["grpBP"][0, g]), peak_nTP=int(pr["grpTP"][0, g]))
                    d1 = describe(net_o, net_n)
                    d1["signflip_p"] = two_sided_p(net_o, fl[:, g]) if net_o != 0 else np.nan
                    ng_rows.append(dict(base, stat="net", **d1))
                    ng_rows.append(dict(base, stat="events_per_snp", **describe(float(eps_o), eps_n)))
                # per-member peak counts (description only)
                col = "direction" if rule == "standard" else "strong_direction"
                sub = calls[(calls["set"] == s) & ((calls["locus"] == sc) if sc != "pooled" else True)]
                for gname in grp["names"]:
                    for mid in NAMED_GROUPS[gname]:
                        mm = sub[sub["motif_id"] == mid]
                        nb, nt = int((mm[col] == "BP-only").sum()), int((mm[col] == "TP-only").sum())
                        name = fam_name_of(hits, mid)
                        ng_rows.append(dict(set=s, rule=rule, scope=sc, group=gname, members=mid,
                                            stat=f"member_count:{name}", peak_nBP=nb, peak_nTP=nt,
                                            observed=nb - nt))
    named = pd.DataFrame(ng_rows)
    named["p_bonferroni"] = np.nan
    for (s, rule, sc, stat), idx in named.groupby(["set", "rule", "scope", "stat"]).groups.items():
        if stat.startswith("member_count"):
            continue
        p = named.loc[idx, "p_two_sided"].to_numpy(dtype=float)
        named.loc[idx, "p_bonferroni"] = np.minimum(1.0, p * len(grp["names"]))
    lead = ["set", "rule", "scope", "group", "stat", "members", "peak_nBP", "peak_nTP", "observed", "null_mean",
            "null_sd", "null_q025", "null_q975", "n_valid_draws", "p_two_sided", "p_bonferroni", "signflip_p"]
    named = named[[c for c in lead if c in named.columns]]
    named.to_csv(out / "named_group_tests.tsv", sep="\t", index=False)

    # ---- calibration and SD guard
    cal_rows, abort_msgs = [], []
    for si, s in enumerate(sets):
        for ri, rule in enumerate(RULES):
            for ci, sc in enumerate(scopes):
                f_null, ev = null_f[(s, rule, sc)]
                p_raw, p_rand = calibration_pvalues(f_null, args.n_calib, _rng(args.seed, 3, si, ri, ci))
                for name, pv in (("ks_uniform_raw", p_raw), ("ks_uniform_randomized", p_rand)):
                    if len(pv) >= 5:
                        ks = stats.kstest(pv, "uniform")
                        cal_rows.append(dict(check=name, set=s, rule=rule, scope=sc, n=len(pv), stat=ks.statistic,
                                             p_value=ks.pvalue, threshold=0.01, passed=bool(ks.pvalue >= 0.01)))
                g = sd_guard(f_null, ev)
                cal_rows.append(dict(check="sd_guard", set=s, rule=rule, scope=sc, n=int((~np.isnan(f_null)).sum()),
                                     stat=g["sd"], threshold=g["threshold"], passed=g["passed"],
                                     detail=f"mean events {g['mean_events']:.2f}"))
                if g["passed"] is False:
                    msg = (f"null SD of f = {g['sd']:.4f} < {SD_GUARD_FACTOR} * 0.5 / sqrt({g['mean_events']:.1f}) = "
                           f"{g['threshold']:.4f} ({s}/{rule}/{sc})")
                    if (s, rule) == PRIMARY:
                        abort_msgs.append(msg)
                    else:
                        warnings.warn("SD guard (sensitivity analysis): " + msg)

    # ---- power (primary analysis, family level)
    power = power_checks(A, sets, loci, fam, args)
    cal_rows += power["rows"]
    calibration = pd.DataFrame(cal_rows)
    calibration.to_csv(out / "calibration.tsv", sep="\t", index=False)
    if abort_msgs:
        (out / "summary.txt").write_text("ABORTED: the matched null failed the SD guard (README.md statistic 6).\n"
                                         + "\n".join(abort_msgs) + "\nSee calibration.tsv and fallbacks.tsv.\n")
        raise CalibrationError("SD guard failed: " + "; ".join(abort_msgs))

    # ---- balance
    bal_rows = []
    for s in sets:
        for sc in scopes:
            pk, nl = A[(s, sc, "peak")], A[(s, sc, "null")]
            for name, fn in balance_covariates():
                o, nv = fn(pk)[0], fn(nl)
                bal_rows.append(dict(set=s, scope=sc, covariate=name, peak=float(o), null_mean=float(np.nanmean(nv)),
                                     null_sd=float(np.nanstd(nv, ddof=1)), peak_units=pk["n_units"],
                                     peak_snps=int(pk["nsnp"][0])))
    balance = pd.DataFrame(bal_rows)
    balance.to_csv(out / "balance.tsv", sep="\t", index=False)

    # ---- plots
    plotted = []
    if not args.no_plots:
        for sc in scopes:
            f_null, _ = null_f[("extended", "standard", sc)]
            row = locus_tests[(locus_tests.set == "extended") & (locus_tests.rule == "standard")
                              & (locus_tests.scope == sc) & (locus_tests.stat == "f")]
            frow_ = locus_tests[(locus_tests.set == "extended") & (locus_tests.rule == "standard")
                                & (locus_tests.scope == sc) & (locus_tests.stat == "f") & (locus_tests.null == "sign_flip")]
            row = row[row.null == "matched"]
            ok = plot_null(out / f"null_fraction_{sc}.png", sc, row, f_null, flips[("extended", "standard", sc)]["f"],
                           float(frow_["p_two_sided"].iloc[0]) if len(frow_) else np.nan)
            if ok:
                plotted.append(f"null_fraction_{sc}.png")

    res = dict(units=units, rows=rows, T=T, A=A, draws=draws, flips=flips, null_f=null_f, fam=fam, grp=grp,
               locus_tests=locus_tests, family_tests=family_tests, named=named, calibration=calibration,
               balance=balance, fallbacks=fallbacks, calls=calls, power=power, loci=loci, sets=sets, strata=strata,
               fam_null=fam_null, runtime=time.time() - t0, plotted=plotted, pol_check=pol_check, n_hits=len(hap))
    write_summary(out / "summary.txt", res, args)
    log(f"done in {time.time() - t0:.1f}s; results in {out}")
    return res


def fam_name_of(hits, motif_id):
    i = np.where(hits["motif_ids"] == motif_id)[0]
    return f"{motif_id} {hits['motif_names'][i[0]]}" if len(i) else f"{motif_id} (absent)"


def balance_covariates():
    def snp_mean(key):
        return lambda a: ratio(a[key], a["sum_k"])
    out = [("k", lambda a: np.asarray(a["sum_k"], float) / max(a["n_units"], 1)),
           ("gc", lambda a: np.asarray(a["sum_gc"], float) / max(a["n_units"], 1)),
           ("cpg", snp_mean("sum_cpg")), (AF_LABEL, snp_mean("sum_af")), ("bp_is_ref", snp_mean("sum_bpref"))]
    for i, c in enumerate(CLASSES):
        out.append((f"class_{c}", (lambda i: lambda a: ratio(a["sum_cls"][:, i], a["sum_k"]))(i)))
    return out


def peak_calls(units, hits, hap, rows, fam, grp, pol):
    """One row per directional (peak unit, motif), standard rule; `strong_direction` marks the strong-change subset."""
    u = hits["unit_idx"]
    sel = units["is_peak"][u]
    slot = hap[sel] - rows["hap_off"][u[sel]]
    df = pd.DataFrame({"unit_idx": u[sel], "slot": slot, "motif_idx": hits["motif_idx"][sel],
                       "mlog10p": hits["mlog10p"][sel]})
    cols = ["locus", "set", "unit_idx", "chr", "k", "pos", "snp_id", "ref", "alt", "bp_alleles", "tp_alleles", "log_p",
            "delta_af", "motif_id", "motif_name", "tf_family", "tf_family_all", "tf_class", "named_group", "direction",
            "strong_direction", "bp_mlog10p", "tp_mlog10p", "_unit_order"]
    if df.empty:
        return pd.DataFrame(columns=cols)
    wide = df.pivot_table(index=["unit_idx", "motif_idx"], columns="slot", values="mlog10p", aggfunc="max")
    wide = wide.reindex(columns=[0, 1]).fillna(0.0).reset_index()
    bpv, tpv = wide[0].to_numpy(), wide[1].to_numpy()
    pb, pt = bpv > PRESENT_MLOG10P, tpv > PRESENT_MLOG10P
    wb, wt = bpv >= SCAN_MLOG10P - MLOG10P_TOL, tpv >= SCAN_MLOG10P - MLOG10P_TOL
    keep = pb != pt
    w = wide[keep].copy()
    w["direction"] = np.where(pb[keep], "BP-only", "TP-only")
    strong = np.where(pb[keep], ~wt[keep], ~wb[keep])
    w["strong_direction"] = np.where(strong, w["direction"], "")
    w["bp_mlog10p"], w["tp_mlog10p"] = bpv[keep], tpv[keep]
    d = units["df"]
    letters = np.array(list(BASES))
    dal = None
    if pol is not None and "delta_af" in pol.columns:
        dal = dict(zip(zip(pol["chr"], pol["pos"]), pol["delta_af"]))
    info = {}
    for uu in np.unique(w["unit_idx"]):
        bp, tp, _ = peak_bases(units, uu)
        pos = d.at[uu, "pos"]
        da = ""
        if dal is not None:
            da = ",".join(str(dal.get((d.at[uu, "chr"], int(x)), "NA")) for x in pos.split(","))
        lp = [float(x) for x in d.at[uu, "log_p"].split(",")]
        info[uu] = dict(locus=d.at[uu, "locus"], set=d.at[uu, "set"], chr=d.at[uu, "chr"], k=int(units["k"][uu]),
                        pos=pos, snp_id=d.at[uu, "snp_id"], ref=d.at[uu, "ref"], alt=d.at[uu, "alt"],
                        bp_alleles=",".join(letters[bp]), tp_alleles=",".join(letters[tp]), log_p=d.at[uu, "log_p"],
                        delta_af=da, _unit_order=-max(lp))
    meta = pd.DataFrame([dict(unit_idx=uu, **v) for uu, v in info.items()])
    w = w.merge(meta, on="unit_idx", how="left")
    mi = w["motif_idx"].to_numpy()
    w["motif_id"] = hits["motif_ids"][mi]
    w["motif_name"] = hits["motif_names"][mi]
    w["tf_family"] = fam["primary"][mi]
    w["tf_family_all"] = fam["all_families"][mi]
    w["tf_class"] = fam["motif_class"][mi]
    gidx = grp["motif_grp"][mi]
    w["named_group"] = np.where(gidx >= 0, np.array(grp["names"] + [""])[gidx], "")
    w["set_order"] = w["set"].map({"extended": 0, "core": 1})
    w = w.sort_values(["set_order", "_unit_order", "unit_idx", "direction", "motif_id"])
    return w[cols].reset_index(drop=True)


def power_checks(A, sets, loci, fam, args):
    """Statistic 7, primary analysis (extended, standard), family level."""
    rows, detail = [], {}
    s, rule = PRIMARY
    for li, L in enumerate(loci):
        pk, nl = A[(s, L, "peak")][rule], A[(s, L, "null")][rule]
        obs = (pk["famBP"] - pk["famTP"])[0].astype(float)
        nul = (nl["famBP"] - nl["famTP"]).astype(float)
        testable = (obs != 0) | (nul != 0).any(0)
        sd = nul.std(0, ddof=1)
        cand = np.where(testable & (sd > 0))[0]
        if len(cand) == 0:
            rows.append(dict(check="power", set=s, rule=rule, scope=L, detail="no testable family", passed=None))
            continue
        order = cand[np.lexsort((np.array(fam["names"])[cand], sd[cand]))]
        target = int(order[len(order) // 2])

        def planted_q(fi, base, null, exclude=None):
            x = base.copy()
            x[fi] += PLANT_N
            if exclude is None:
                p = two_sided_p(x, null)
            else:   # leave-one-out for a pseudo-peak taken from the null
                nv = null.shape[0] - 1
                lo = (null <= x).sum(0) - (null[exclude] <= x)
                hi = (null >= x).sum(0) - (null[exclude] >= x)
                p = np.minimum(1.0, 2 * np.minimum((lo + 1) / (nv + 1), (hi + 1) / (nv + 1)))
            p = np.where(testable, p, np.nan)
            return bh_fdr(p)[fi]

        q_obs = planted_q(target, obs, nul)
        sweep = np.array([planted_q(fi, obs, nul) for fi in cand])
        npw = min(args.n_power, nul.shape[0])
        qs = np.array([planted_q(target, nul[d], nul, exclude=d) for d in range(npw)])
        name = fam["names"][target]
        detail[L] = dict(family=name, null_sd=float(sd[target]), q_observed=float(q_obs),
                         power=float((qs < 0.05).mean()), median_q=float(np.median(qs)),
                         sweep_pass=float((sweep < 0.05).mean()), n_tested=len(cand),
                         max_sd_pass=float(sd[cand][sweep < 0.05].max()) if (sweep < 0.05).any() else np.nan)
        rows.append(dict(check="power_observed_plus20", set=s, rule=rule, scope=L, n=len(cand), stat=PLANT_N,
                         p_value=q_obs, threshold=0.05, passed=bool(q_obs < 0.05),
                         detail=f"family {name} (median null SD {sd[target]:.2f}); value = q"))
        rows.append(dict(check="power_pseudo_peaks", set=s, rule=rule, scope=L, n=npw, stat=float((qs < 0.05).mean()),
                         p_value=float(np.median(qs)), threshold=0.05, passed=bool(np.median(qs) < 0.05),
                         detail=f"family {name}; stat = share of pseudo-peaks with q < 0.05, value = median q"))
        rows.append(dict(check="power_sweep_all_families", set=s, rule=rule, scope=L, n=len(cand),
                         stat=float((sweep < 0.05).mean()), threshold=0.05,
                         detail=f"share of testable families reaching q < 0.05 with +{PLANT_N}; largest null SD that "
                                f"passes {detail[L]['max_sd_pass']:.2f}"))
    return dict(rows=rows, detail=detail)


def plot_null(path, scope, row, f_null, f_flip, p_flip=np.nan):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        warnings.warn("matplotlib missing; plots skipped")
        return False
    if row.empty or np.isnan(row["observed"].iloc[0]):
        return False
    obs = float(row["observed"].iloc[0])
    p = float(row["p_two_sided"].iloc[0])
    ink, null_col, grid = "#1f2328", "#7a9cc6", "#d0d7de"
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.2), sharex=True)
    for ax, vals, title in ((axes[0], f_null, "Matched random SNPs"), (axes[1], f_flip, "Sign-flip of peak units")):
        v = np.asarray(vals, float)
        v = v[~np.isnan(v)]
        ax.hist(v, bins=40, color=null_col, edgecolor="white", linewidth=0.5)
        ax.axvline(obs, color=ink, linewidth=2)
        ax.axvline(0.5, color=grid, linewidth=1, linestyle="--", zorder=0)
        ax.set_title(title, fontsize=10, color=ink, loc="left")
        ax.set_xlabel("BP-only fraction f = nBP / (nBP + nTP)", fontsize=9, color=ink)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(labelsize=8, colors=ink)
    axes[0].set_ylabel("null draws", fontsize=9, color=ink)
    for ax, pv in ((axes[0], p), (axes[1], p_flip)):
        ax.annotate(f"peak f = {obs:.3f}\np = {pv:.3g}", xy=(obs, ax.get_ylim()[1] * 0.92), xytext=(6, 0),
                    textcoords="offset points", fontsize=8, color=ink, va="top")
    fig.suptitle(f"{scope}: peak BP-only fraction vs nulls (set extended)", fontsize=10, color=ink, x=0.01,
                 ha="left")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return True


def _short_list(s, n=10):
    xs = s.split(";")
    return ";".join(xs[:n]) + (f"; +{len(xs) - n} more (locus_tests.tsv)" if len(xs) > n else "")


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "NA"
    return f"{x:.{nd}g}" if abs(x) < 1e-3 and x != 0 else f"{x:.{nd}f}"


def write_summary(path, res, args):
    lt, ft, nm, cal, bal, fb = (res["locus_tests"], res["family_tests"], res["named"], res["calibration"],
                                res["balance"], res["fallbacks"])
    units, strata, loci, sets = res["units"], res["strata"], res["loci"], res["sets"]
    L = []
    w = L.append
    w("TFBS matched-null test (code/15_TFBS_Scan/README.md)")
    w(f"Run {time.strftime('%Y-%m-%d %H:%M')}; seed {args.seed}; {args.n_draws} draws per locus and set; "
      f"{args.n_calib} calibration pseudo-peaks; min stratum {args.min_stratum}; runtime {res['runtime']:.0f} s")
    w(f"Inputs: units {args.units}; hits {args.hits}; families {args.families}; polarized {args.polarized or 'not given'}")
    pk = units["is_peak"]
    k = units["k"]
    w(f"Units: {int(pk.sum())} peak, {strata['n_pool1']} pool singletons, {strata['n_pool2']} pool chains "
      f"(k: {', '.join(f'{i}={c}' for i, c in enumerate(strata['pool2_k']) if c)}); KMAX {strata['kmax']}; "
      f"{res['n_hits']} hits; {len(res['fam']['names'])} families")
    if res["pol_check"]:
        w(f"bp_mask agrees with peak_snps_polarized.tsv at {res['pol_check'][0]} peak SNPs "
          f"({res['pol_check'][1]} not found there)")
    w(f"Pool GC quintile edges (k=1): {np.round(strata['gc1_edges'], 4).tolist()}; {AF_LABEL} quintile edges: "
      f"{np.round(strata['af_edges'], 4).tolist()}")
    if strata["gc2_edges"] is not None:
        w(f"Pool GC quintile edges (chains): {np.round(strata['gc2_edges'], 4).tolist()}")
    w("")
    w("Peak units per locus (set: units / SNPs / max k):")
    for s in sets:
        for Lc in loci:
            sel = pk & (units["set"] == s) & (units["locus"] == Lc)
            w(f"  {s:8s} {Lc:8s} {int(sel.sum()):4d} units / {int(k[sel].sum()):4d} SNPs / max k "
              f"{int(k[sel].max()) if sel.any() else 0}")
    w("Stratum fallbacks (fallbacks.tsv):")
    for s in sets:
        sub = fb[fb["set"] == s]
        cnt = sub.groupby("level").size().to_dict()
        w(f"  {s}: " + ", ".join(f"{a} {b}" for a, b in cnt.items()) + f"; below min {int(sub['below_min'].sum())}; "
          f"median pool size {int(sub['pool_size'].median()) if len(sub) else 0}; groups drawn with replacement "
          f"{int(sub['with_replacement'].sum()) if len(sub) else 0}")
    w("")

    def frow(s, rule, sc, stat="f", null="matched"):
        r = lt[(lt.set == s) & (lt.rule == rule) & (lt.scope == sc) & (lt.stat == stat) & (lt.null == null)]
        return r.iloc[0] if len(r) else None

    w("1. PRIMARY: pooled BP-only fraction f, set extended, standard rule (motif-level calls)")
    r = frow("extended", "standard", "pooled")
    w(f"   peak f = {_fmt(r.observed)} (nBP {int(r.peak_nBP)}, nTP {int(r.peak_nTP)}); matched null mean {_fmt(r.null_mean)} "
      f"(SD {_fmt(r.null_sd)}, 95% {_fmt(r.null_q025)}-{_fmt(r.null_q975)}); two-sided p = {_fmt(r.p_two_sided)}")
    w("2. Per locus (extended, standard):")
    for sc in loci:
        r = frow("extended", "standard", sc)
        e = frow("extended", "standard", sc, "events_per_snp")
        w(f"   {sc:8s} f = {_fmt(r.observed)} (nBP {int(r.peak_nBP)}, nTP {int(r.peak_nTP)}) vs null {_fmt(r.null_mean)} "
          f"[{_fmt(r.null_q025)}, {_fmt(r.null_q975)}], p = {_fmt(r.p_two_sided)}; events/SNP {_fmt(e.observed)} vs "
          f"{_fmt(e.null_mean)} [{_fmt(e.null_q025)}, {_fmt(e.null_q975)}], p = {_fmt(e.p_two_sided)}")
    r = frow("extended", "standard", "pooled", "events_per_snp")
    w(f"   pooled   events/SNP {_fmt(r.observed)} vs {_fmt(r.null_mean)}, p = {_fmt(r.p_two_sided)}")
    r = frow("extended", "standard", "pooled", "f_family")
    w(f"   (secondary) family-level pooled f = {_fmt(r.observed)} vs {_fmt(r.null_mean)}, p = {_fmt(r.p_two_sided)}")
    r = frow("extended", "standard", "pooled", "f_equal_weight")
    w(f"   (secondary) equal-weight pooled f (mean of locus f) = {_fmt(r.observed)} vs {_fmt(r.null_mean)} "
      f"[{_fmt(r.null_q025)}, {_fmt(r.null_q975)}], p = {_fmt(r.p_two_sided)}")
    w("3. Families (extended, standard; BH within each scope):")
    sub = ft[(ft.set == "extended") & (ft.rule == "standard")] if len(ft) else ft
    for sc in loci + ["pooled"]:
        s2 = sub[sub.scope == sc] if len(sub) else sub
        hits = s2[s2.q_bh < 0.05] if len(s2) else s2
        w(f"   {sc:8s} {len(s2)} testable families; q < 0.05: {len(hits)}"
          + ("" if hits.empty else "; " + "; ".join(f"{x.tf_family} net {x.peak_net:+d} (null {x.null_mean_net:+.2f}, "
                                                       f"q {_fmt(x.q_bh)})" for x in hits.head(10).itertuples())))
    w("4. Recurrence (families with nBP - nTP > 0 at >= 2 loci; Unknown excluded):")
    for s in sets:
        for rule in RULES:
            r = frow(s, rule, "all_loci", "recurrence", "matched_joint")
            if r is not None:
                w(f"   {s:8s} {rule:8s} observed {int(r.observed)} vs null mean {_fmt(r.null_mean)} "
                  f"[{_fmt(r.null_q025)}, {_fmt(r.null_q975)}], p = {_fmt(r.p_two_sided)}"
                  + (f" ({_short_list(r.detail)})" if isinstance(r.detail, str) and r.detail else ""))
    w("4b. Named groups (Bonferroni x2 within scope; outside the family BH):")
    for s in sets:
        for rule in RULES:
            for sc in loci + ["pooled"]:
                for gname in res["grp"]["names"]:
                    x = nm[(nm.set == s) & (nm.rule == rule) & (nm.scope == sc) & (nm.group == gname)]
                    a = x[x.stat == "net"]
                    b = x[x.stat == "events_per_snp"]
                    if a.empty:
                        continue
                    a, b = a.iloc[0], b.iloc[0]
                    if (s, rule) != PRIMARY and sc != "pooled":
                        continue
                    mem = x[x.stat.str.startswith("member_count")]
                    memtxt = ", ".join(f"{m.stat.split(':', 1)[1]} {m.peak_nBP}/{m.peak_nTP}" for m in mem.itertuples())
                    w(f"   {s:8s} {rule:8s} {sc:8s} {gname:4s} net {int(a.observed):+d} (BP {a.peak_nBP}, TP "
                      f"{a.peak_nTP}) vs null {_fmt(a.null_mean)} [{_fmt(a.null_q025)}, {_fmt(a.null_q975)}], p = "
                      f"{_fmt(a.p_two_sided)}, Bonferroni {_fmt(a.p_bonferroni)}; events/SNP {_fmt(b.observed)} vs "
                      f"{_fmt(b.null_mean)}, p = {_fmt(b.p_two_sided)}, Bonferroni {_fmt(b.p_bonferroni)}"
                      + (f"; members BP/TP: {memtxt}" if (s, rule) == PRIMARY else ""))
    w("5. Sensitivity (pooled f):")
    r = frow("extended", "standard", "pooled", "f", "sign_flip")
    w(f"   (a) sign-flip null: mean {_fmt(r.null_mean)} [{_fmt(r.null_q025)}, {_fmt(r.null_q975)}], "
      f"p = {_fmt(r.p_two_sided)}")
    for sc in loci:
        r = frow("extended", "standard", sc, "f", "sign_flip")
        w(f"       {sc:8s} sign-flip p = {_fmt(r.p_two_sided)}")
    r = frow("extended", "strong", "pooled")
    w(f"   (b) strong-change rule: f = {_fmt(r.observed)} (nBP {int(r.peak_nBP)}, nTP {int(r.peak_nTP)}) vs "
      f"{_fmt(r.null_mean)}, p = {_fmt(r.p_two_sided)}")
    if "core" in sets:
        r = frow("core", "standard", "pooled")
        w(f"   (c) set core: f = {_fmt(r.observed)} (nBP {int(r.peak_nBP)}, nTP {int(r.peak_nTP)}) vs {_fmt(r.null_mean)}, "
          f"p = {_fmt(r.p_two_sided)}")
        for sc in loci:
            r = frow("core", "standard", sc)
            w(f"       {sc:8s} core f = {_fmt(r.observed)} vs {_fmt(r.null_mean)}, p = {_fmt(r.p_two_sided)}")
    w("6. Calibration (extended, standard; 500 pseudo-peaks, KS vs uniform) and SD guard:")
    for sc in loci + ["pooled"]:
        c = cal[(cal.set == "extended") & (cal.rule == "standard") & (cal.scope == sc)]
        ks1 = c[c.check == "ks_uniform_raw"]
        ks2 = c[c.check == "ks_uniform_randomized"]
        g = c[c.check == "sd_guard"]
        txt = f"   {sc:8s}"
        if len(ks1):
            txt += f" KS D = {_fmt(ks1.stat.iloc[0])} (p {_fmt(ks1.p_value.iloc[0])}); tie-randomised D = " \
                   f"{_fmt(ks2.stat.iloc[0])} (p {_fmt(ks2.p_value.iloc[0])});"
        if len(g):
            txt += f" null SD {_fmt(g.stat.iloc[0])} vs guard {_fmt(g.threshold.iloc[0])} " \
                   f"({'pass' if g.passed.iloc[0] else 'FAIL' if g.passed.iloc[0] is False else 'n/a'})"
        w(txt)
    w(f"7. Power (+{PLANT_N} BP-only calls on one family, extended, standard):")
    for Lc, d in res["power"]["detail"].items():
        w(f"   {Lc:8s} {d['family']} (median null SD {d['null_sd']:.2f}): q = {_fmt(d['q_observed'])} "
          f"({'pass' if d['q_observed'] < 0.05 else 'FAIL'}); pseudo-peak power {d['power']:.2f} (median q "
          f"{_fmt(d['median_q'])}); {d['sweep_pass']:.2f} of {d['n_tested']} families reach q < 0.05")
    w("8. Balance (extended; peak vs matched-null mean):")
    for sc in loci + ["pooled"]:
        b = bal[(bal.set == "extended") & (bal.scope == sc)]
        w(f"   {sc:8s} " + "; ".join(f"{x.covariate} {x.peak:.3f}/{x.null_mean:.3f}" for x in b.itertuples()))
    w("")
    w("Outputs: locus_tests.tsv, family_tests.tsv, named_group_tests.tsv, balance.tsv, calibration.tsv, fallbacks.tsv, "
      + ", ".join(f"per_unit_calls_{x}.tsv" for x in loci) + (", " + ", ".join(res["plotted"]) if res["plotted"] else ""))
    w("Choices that README.md leaves open: see the docstring of matched_null.py.")
    Path(path).write_text("\n".join(L) + "\n")


def main(argv=None):
    args = parse_args(argv)
    try:
        run(args)
    except CalibrationError as e:
        print(f"ABORT: {e}", file=sys.stderr)
        print("The matched null is narrower than a binomial; check units/hits before citing any result.",
              file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
