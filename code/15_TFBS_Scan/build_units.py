#!/usr/bin/env python3
"""
Stage A2 of the TFBS scan (specification: code/15_TFBS_Scan/README.md).

Build the scanning units and their haplotype sequences:

  Peak units  non-CDS peak SNPs of one locus (from peak_snps_polarized.tsv),
              sorted by position; consecutive SNPs with gap <= F share a unit.
              set `extended` clusters all non-CDS peak SNPs; set `core`
              clusters only the core (log_p > 5) non-CDS SNPs.
  Pool units  GEMMA SNPs with af > 0.15 (the GWAS filter), single-base ACGT alleles,
              REF (allele1) equal to the FASTA base, outside extended peaks
              +-100 kb (peaks.bed, numeric chr -> 'chr'+n) and outside CDS.
              k = 1: a random sample of pool SNPs (default 300,000).
              k >= 2: random chains of k consecutive pool SNPs with every gap
              <= F, up to 20,000 per k for k = 2..KMAX, KMAX = min(4, max peak k).
              Units whose span +-F leaves the chromosome or contains N are
              removed from the candidate set before sampling, so the sampled
              counts are exact. Chains need not be maximal: other pool SNPs
              within F stay reference, as non-peak SNPs do at the peak.

Outputs (in --out_dir):
  units.tsv  one row per unit (columns as in README.md; peak units first)
  units.fa   peak units: masks bp_mask and bp_mask ^ (2^k - 1); pool units:
             all 2^k masks. Header >U{unit_idx:07d}_M{mask}; bit j = 1 means
             ALT at SNP j. Each sequence is ref[pos_0 - F .. pos_last + F].

Features: gc = GC fraction of the reference over [pos_0 - 50, pos_last + 50]
(N excluded); cpg per SNP = 1 if either allele is C followed by reference G or
G preceded by reference C; gap_bin = minimum gap binned 1-4 / 5-12 / 13-32
(NA for k = 1); af_ref = GEMMA `af` per SNP (see the note in
polarize_peak_snps.py: this `af` tracks the ALT frequency; copied unchanged).

After writing, both files are read back and every unit and sequence is checked
against the reference through an independent random-access FASTA path
(unit_idx 0..N-1, coordinates/offsets, mask-0 = reference substring, REF base
at every snp_off, Hamming distance = popcount(mask) at exactly the snp_off
positions with the right alleles, ACGT only (+N only where the reference has N
in peak units), gc and cpg recomputed). Any failure aborts with an error.

Usage:
  python3 code/15_TFBS_Scan/build_units.py \\
    [--polarized output_data/15_TFBS_Scan/peak_snps_polarized.tsv] \\
    [--assoc output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/PKINGS_ALL_WOB_EXCLUDED.lmm2.assoc.txt] \\
    [--peaks_bed input_data/10_Introgression/peaks.bed] [--ref $TFBS_REF] [--gff $TFBS_GFF] \\
    [--out_dir output_data/15_TFBS_Scan] [--n_singletons 300000] [--n_chains 20000] \\
    [--kmax_cap 4] [--seed 20261013]
"""

import argparse
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from polarize_peak_snps import (  # noqa: E402
    ROOT, IndexedFasta, env_default, in_intervals_vec, read_cds_intervals,
)

W = 33            # widest JASPAR motif
F = W - 1         # flank on each side; also the maximum within-unit gap
GC_FLANK = 50     # gc window: +-50 bp around the unit span
LOCI = ["sptbn4b", "cntnap5", "itih3"]
ACGT = set("ACGT")

TSV_COLUMNS = [
    "unit_idx", "source", "locus", "set", "chr", "k", "pos", "ref", "alt", "snp_off", "seq_start",
    "seq_len", "gap_bin", "gc", "cpg", "af_ref", "snp_id", "log_p", "bp_mask",
]


class Unit:
    """One scanning unit: k SNPs sorted by position plus derived fields."""
    __slots__ = ("idx", "source", "locus", "set", "chr", "pos", "ref", "alt", "af", "snp_id",
                 "log_p", "bp_mask", "gc", "cpg", "seq0")

    def __init__(self, source, locus, set_, chrom, pos, ref, alt, af, snp_id, log_p=None, bp_mask=None):
        self.idx = -1
        self.source, self.locus, self.set, self.chr = source, locus, set_, chrom
        self.pos, self.ref, self.alt, self.af, self.snp_id = pos, ref, alt, af, snp_id
        self.log_p, self.bp_mask = log_p, bp_mask
        self.gc = self.cpg = self.seq0 = None

    @property
    def k(self):
        return len(self.pos)

    @property
    def seq_start(self):
        return self.pos[0] - F

    @property
    def seq_len(self):
        return self.pos[-1] - self.pos[0] + 2 * F + 1

    @property
    def snp_off(self):
        return [p - self.seq_start + 1 for p in self.pos]


def gap_bin(pos: list) -> str:
    if len(pos) < 2:
        return "NA"
    g = min(b - a for a, b in zip(pos, pos[1:]))
    if not 1 <= g <= F:
        raise ValueError(f"min gap {g} outside 1..{F}")
    return "1-4" if g <= 4 else ("5-12" if g <= 12 else "13-32")


def gc_fraction(s: str) -> float:
    n_acgt = len(s) - s.count("N")
    return (s.count("G") + s.count("C")) / n_acgt if n_acgt else float("nan")


def cpg_flags(seq0: str, offs: list, refs: list, alts: list) -> list:
    """1 if either allele forms a CpG with the reference neighbour (offs 1-based in seq0)."""
    out = []
    for o, r, a in zip(offs, refs, alts):
        prev, nxt = seq0[o - 2], seq0[o]
        out.append(int(any((x == "C" and nxt == "G") or (x == "G" and prev == "C") for x in (r, a))))
    return out


def annotate(unit: Unit, get) -> None:
    """Fill seq0 (mask-0 sequence), gc and cpg; `get(start, end)` returns 1-based inclusive ref."""
    unit.seq0 = get(unit.seq_start, unit.seq_start + unit.seq_len - 1)
    lo, hi = max(1, unit.pos[0] - GC_FLANK), min(get.length, unit.pos[-1] + GC_FLANK)
    unit.gc = gc_fraction(get(lo, hi))
    unit.cpg = cpg_flags(unit.seq0, unit.snp_off, unit.ref, unit.alt)


def apply_mask(seq0: str, offs: list, alts: list, mask: int) -> str:
    s = list(seq0)
    for j, (o, a) in enumerate(zip(offs, alts)):
        if mask >> j & 1:
            s[o - 1] = a
    return "".join(s)


# ── Peak units ────────────────────────────────────────────────────────────────
def clusters(pos: list) -> list:
    """Index groups of sorted positions; consecutive positions with gap <= F share a group."""
    groups = []
    for i, p in enumerate(pos):
        if groups and p - pos[i - 1] <= F:
            groups[-1].append(i)
        else:
            groups.append([i])
    return groups


def build_peak_units(snps: pd.DataFrame) -> tuple:
    if not snps["ref_ok"].all():
        sys.exit(f"ERROR: {int((~snps['ref_ok']).sum())} peak SNPs have REF != FASTA; fix stage A1 first")
    bad = ~(snps["ref"].isin(ACGT) & snps["alt"].isin(ACGT) & (snps["ref"] != snps["alt"]))
    if bad.any():
        sys.exit(f"ERROR: non-ACGT or REF == ALT peak alleles:\n{snps[bad]}")
    units, cds_dropped = [], {}
    for set_ in ("extended", "core"):
        for locus in LOCI:
            t = snps[snps["locus"] == locus].sort_values("pos")
            cds_dropped[locus] = int(t["in_cds"].sum())
            t = t[~t["in_cds"]]
            if set_ == "core":
                t = t[t["core"]]
            pos = t["pos"].tolist()
            if len(set(pos)) != len(pos):
                sys.exit(f"ERROR: duplicate peak positions at {locus}")
            for g in clusters(pos):
                r = t.iloc[g]
                bp_mask = sum(1 << j for j, is_ref in enumerate(r["bp_is_ref"]) if not is_ref)
                units.append(Unit("peak", locus, set_, r["chr"].iloc[0], r["pos"].tolist(),
                                  r["ref"].tolist(), r["alt"].tolist(), r["af_ref"].tolist(),
                                  r["snp_id"].tolist(), r["log_p"].tolist(), bp_mask))
    return units, cds_dropped


# ── Pool ──────────────────────────────────────────────────────────────────────
def load_pool(assoc: str, af_min: float, fasta: IndexedFasta) -> tuple:
    """Filtered GEMMA SNPs (chr, ps, allele1, allele0, af), sorted in FASTA order."""
    counts = {}
    df = pd.read_csv(assoc, sep="\t", usecols=["chr", "ps", "allele1", "allele0", "af"],
                     dtype={"chr": "category", "ps": np.int64, "allele1": "category",
                            "allele0": "category", "af": np.float64})
    counts["rows in assoc"] = len(df)
    df = df[df["af"] > af_min]
    counts[f"af > {af_min}"] = len(df)
    df = df[df["allele1"].isin(ACGT) & df["allele0"].isin(ACGT) & (df["allele1"].astype(str)
                                                                   != df["allele0"].astype(str))]
    counts["single-base ACGT alleles"] = len(df)
    dup = df.duplicated(["chr", "ps"], keep=False)
    df = df[~dup]
    counts["no duplicate positions"] = len(df)
    unknown = set(df["chr"].astype(str).unique()) - set(fasta.index)
    if unknown:
        sys.exit(f"ERROR: assoc chromosomes missing from FASTA index: {sorted(unknown)[:10]}")
    order = {c: i for i, c in enumerate(fasta.order)}
    df = df.assign(chr=df["chr"].astype(str))
    df["chr_code"] = df["chr"].map(order).astype(np.int32)
    df = df.sort_values(["chr_code", "ps"], kind="stable").reset_index(drop=True)
    return df, counts


def read_peaks_bed(path: str, buffer: int) -> dict:
    """{chrom: [(lo, hi), ...]} exclusion intervals (1-based inclusive), peaks +- buffer."""
    bed = pd.read_csv(path, sep="\t", header=None)
    out = {}
    for c, s, e in zip(bed[0], bed[1], bed[2]):
        chrom = str(c) if str(c).startswith("chr") else f"chr{c}"
        out.setdefault(chrom, []).append((int(s) - buffer, int(e) + buffer))
    return out


def filter_pool(df: pd.DataFrame, peaks: dict, cds: dict, fasta: IndexedFasta, counts: dict) -> tuple:
    """Drop peak-proximal, CDS and REF-mismatch SNPs; return (pool, N runs per chrom)."""
    keep = np.ones(len(df), dtype=bool)
    pos_all = df["ps"].to_numpy()
    chr_all = df["chr"].to_numpy()
    for chrom, ivs in peaks.items():
        m = chr_all == chrom
        for lo, hi in ivs:
            keep[m & (pos_all >= lo) & (pos_all <= hi)] = False
    counts["outside extended peaks +-buffer"] = int(keep.sum())
    for chrom, (starts, ends) in cds.items():
        m = np.flatnonzero(chr_all == chrom)
        if len(m):
            keep[m[in_intervals_vec(pos_all[m], starts, ends)]] = False
    counts["outside CDS"] = int(keep.sum())

    n_runs, mismatch = {}, 0
    a1 = np.frombuffer("".join(df["allele1"].astype(str)).encode("ascii"), dtype=np.uint8)
    for chrom in fasta.order:
        m = np.flatnonzero((chr_all == chrom) & keep)
        arr = None
        if len(m):
            arr = fasta.chrom_array(chrom)
            ref_base = arr[pos_all[m] - 1]
            bad = ref_base != a1[m]
            mismatch += int(bad.sum())
            keep[m[bad]] = False
        if arr is None:
            continue
        non_acgt = ~np.isin(arr, np.frombuffer(b"ACGT", dtype=np.uint8))
        idx = np.flatnonzero(non_acgt) + 1  # 1-based positions
        if len(idx):
            brk = np.flatnonzero(np.diff(idx) > 1)
            starts = np.concatenate([[idx[0]], idx[brk + 1]])
            ends = np.concatenate([idx[brk], [idx[-1]]])
        else:
            starts = ends = np.zeros(0, dtype=np.int64)
        n_runs[chrom] = (starts, ends)
    counts["REF (allele1) == FASTA base"] = int(keep.sum())
    counts["  (REF mismatches dropped)"] = mismatch
    return df[keep].reset_index(drop=True), n_runs


def span_ok(pool: pd.DataFrame, i: np.ndarray, j: np.ndarray, n_runs: dict, fasta: IndexedFasta) -> np.ndarray:
    """span [pos_i - F, pos_j + F] inside the chromosome and free of non-ACGT bases."""
    chrom = pool["chr"].to_numpy()[i]
    lo = pool["ps"].to_numpy()[i] - F
    hi = pool["ps"].to_numpy()[j] + F
    ok = lo >= 1
    for c in np.unique(chrom):
        m = np.flatnonzero(chrom == c)
        ok[m] &= hi[m] <= fasta.length(c)
        starts, ends = n_runs.get(c, (np.zeros(0, np.int64), np.zeros(0, np.int64)))
        if len(starts):
            r = np.searchsorted(ends, lo[m], side="left")  # first N run ending at/after lo
            has_n = (r < len(starts)) & (starts[np.minimum(r, len(starts) - 1)] <= hi[m])
            ok[m] &= ~has_n
    return ok


def sample_pool_units(pool: pd.DataFrame, n_runs: dict, fasta: IndexedFasta, kmax: int,
                      n_singletons: int, n_chains: int, rng: np.random.Generator) -> tuple:
    n = len(pool)
    pos, code = pool["ps"].to_numpy(), pool["chr_code"].to_numpy()
    gap_ok = (code[1:] == code[:-1]) & (pos[1:] - pos[:-1] <= F)
    bad_cum = np.concatenate([[0], np.cumsum(~gap_ok)])
    starts_by_k, stats = {}, []
    for k in range(1, kmax + 1):
        i = np.arange(n - k + 1)
        chain = bad_cum[i + k - 1] - bad_cum[i] == 0
        n_chain = int(chain.sum())
        i = i[chain]
        i = i[span_ok(pool, i, i + k - 1, n_runs, fasta)]
        target = n_singletons if k == 1 else n_chains
        take = min(target, len(i))
        starts_by_k[k] = np.sort(rng.choice(i, size=take, replace=False))
        stats.append((k, n_chain, len(i), take))
    units = []
    chr_s, ref_s, alt_s = pool["chr"].to_numpy(), pool["allele1"].astype(str).to_numpy(), \
        pool["allele0"].astype(str).to_numpy()
    af = pool["af"].to_numpy()
    for k in range(1, kmax + 1):
        for s in starts_by_k[k]:
            sl = slice(s, s + k)
            c = chr_s[s]
            p = pos[sl].tolist()
            units.append(Unit("pool", "pool", "pool", c, p, ref_s[sl].tolist(), alt_s[sl].tolist(),
                              af[sl].tolist(), [f"{c}_{x}" for x in p]))
    return units, stats


# ── Output ────────────────────────────────────────────────────────────────────
def unit_masks(u: Unit) -> list:
    full = (1 << u.k) - 1
    return [u.bp_mask, u.bp_mask ^ full] if u.source == "peak" else list(range(1 << u.k))


def fmt_list(xs) -> str:
    return ",".join(str(x) for x in xs)


def write_outputs(units: list, tsv_path: Path, fa_path: Path) -> None:
    with open(tsv_path, "w") as tsv, open(fa_path, "w") as fa:
        tsv.write("\t".join(TSV_COLUMNS) + "\n")
        for u in units:
            row = [u.idx, u.source, u.locus, u.set, u.chr, u.k, fmt_list(u.pos), fmt_list(u.ref),
                   fmt_list(u.alt), fmt_list(u.snp_off), u.seq_start, u.seq_len, gap_bin(u.pos),
                   f"{u.gc:.4f}", fmt_list(u.cpg), fmt_list(u.af), fmt_list(u.snp_id),
                   fmt_list(u.log_p) if u.log_p is not None else "NA",
                   u.bp_mask if u.bp_mask is not None else "NA"]
            tsv.write("\t".join(str(x) for x in row) + "\n")
            offs = u.snp_off
            for m in unit_masks(u):
                fa.write(f">U{u.idx:07d}_M{m}\n{apply_mask(u.seq0, offs, u.alt, m)}\n")


def check(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def validate(tsv_path: Path, fa_path: Path, fasta: IndexedFasta) -> tuple:
    """Re-read units.tsv / units.fa and check every unit and sequence; returns (n_seq, total_bp)."""
    t = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
    check(list(t.columns) == TSV_COLUMNS, "units.tsv columns differ from the specification")
    idx = t["unit_idx"].astype(int).to_numpy()
    check(np.array_equal(idx, np.arange(len(t))), "unit_idx is not 0..N-1 in order")
    src = t["source"].to_numpy()
    n_peak = int((src == "peak").sum())
    check((src[:n_peak] == "peak").all() and (src[n_peak:] == "pool").all(), "peak units are not first")

    seqs = {}
    with open(fa_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                header = line[1:]
                check(header not in seqs, f"duplicate FASTA header {header}")
                seqs[header] = ""
            else:
                seqs[header] += line
    n_bp = sum(len(s) for s in seqs.values())

    expected = 0
    for r in t.itertuples(index=False):
        u = int(r.unit_idx)
        k = int(r.k)
        pos = [int(x) for x in r.pos.split(",")]
        ref, alt = r.ref.split(","), r.alt.split(",")
        offs = [int(x) for x in r.snp_off.split(",")]
        start, slen = int(r.seq_start), int(r.seq_len)
        tag = f"U{u:07d}"
        check(len(pos) == len(ref) == len(alt) == len(offs) == k
              == len(r.cpg.split(",")) == len(r.af_ref.split(",")) == len(r.snp_id.split(",")),
              f"{tag}: list lengths differ from k")
        check(all(b - a >= 1 and b - a <= F for a, b in zip(pos, pos[1:])), f"{tag}: gap outside 1..F")
        check(start == pos[0] - F and offs[0] == F + 1, f"{tag}: seq_start / first snp_off wrong")
        check(slen == pos[-1] - pos[0] + 2 * F + 1, f"{tag}: seq_len wrong")
        check(offs == [p - start + 1 for p in pos], f"{tag}: snp_off wrong")
        check(r.gap_bin == gap_bin(pos), f"{tag}: gap_bin wrong")
        check(all(a in ACGT and b in ACGT and a != b for a, b in zip(ref, alt)), f"{tag}: bad alleles")
        ref0 = fasta.fetch(r.chr, start, start + slen - 1)
        check(all(ref0[o - 1] == a for o, a in zip(offs, ref)), f"{tag}: REF != reference at snp_off")
        check([str(x) for x in cpg_flags(ref0, offs, ref, alt)] == r.cpg.split(","), f"{tag}: cpg wrong")
        lo, hi = max(1, pos[0] - GC_FLANK), min(fasta.length(r.chr), pos[-1] + GC_FLANK)
        check(abs(gc_fraction(fasta.fetch(r.chr, lo, hi)) - float(r.gc)) < 1e-4, f"{tag}: gc wrong")
        full = (1 << k) - 1
        if r.source == "peak":
            check(r.locus in LOCI and r.set in ("extended", "core") and r.log_p != "NA", f"{tag}: labels")
            bp = int(r.bp_mask)
            check(0 <= bp <= full, f"{tag}: bp_mask out of range")
            masks = [bp, bp ^ full]
            allowed = ACGT | {"N"}
        else:
            check(r.locus == "pool" and r.set == "pool" and r.log_p == "NA" and r.bp_mask == "NA",
                  f"{tag}: pool labels")
            check(set(ref0) <= ACGT, f"{tag}: pool unit span +-F contains non-ACGT")
            masks = list(range(full + 1))
            allowed = ACGT
        for m in masks:
            h = f"{tag}_M{m}"
            check(h in seqs, f"missing sequence {h}")
            s = seqs[h]
            expected += 1
            check(len(s) == slen, f"{h}: length {len(s)} != seq_len {slen}")
            check(set(s) <= allowed, f"{h}: unexpected characters")
            diff = [i for i, (a, b) in enumerate(zip(s, ref0)) if a != b]
            want = [offs[j] - 1 for j in range(k) if m >> j & 1]
            check(diff == want and len(diff) == bin(m).count("1"),
                  f"{h}: differs from mask 0 at {diff}, expected {want}")
            check(all(s[offs[j] - 1] == (alt[j] if m >> j & 1 else ref[j]) for j in range(k)),
                  f"{h}: wrong allele at a snp_off")
            if m == 0:
                check(s == ref0, f"{h}: mask-0 sequence != reference substring")
            check(all(c == "N" and ref0[i] == "N" for i, c in enumerate(s) if c not in ACGT),
                  f"{h}: non-ACGT base that is not a reference N")
    check(expected == len(seqs), f"units.fa has {len(seqs) - expected} unexpected sequences")
    return len(seqs), n_bp


# ── Main ──────────────────────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build peak and pool scanning units (TFBS stage A2).")
    p.add_argument("--polarized", default=str(ROOT / "output_data/15_TFBS_Scan/peak_snps_polarized.tsv"),
                   help="Stage A1 output")
    p.add_argument("--assoc", default=str(ROOT / "output_data/06_Association/PKINGS_ALL_WOB_EXCLUDED/"
                                                 "PKINGS_ALL_WOB_EXCLUDED.lmm2.assoc.txt"),
                   help="GEMMA lmm2 association file (pool)")
    p.add_argument("--peaks_bed", default=str(ROOT / "input_data/10_Introgression/peaks.bed"),
                   help="Extended peaks (numeric chromosomes)")
    p.add_argument("--ref", default=env_default(
        "TFBS_REF", "input_data/00_Reference_Genome/"
        "jordan-mic4273-mb-hirise-20xx3__06-13-2024__final_assembly.fasta.renamed.fa"),
        help="Reference FASTA with .fai (default: $TFBS_REF)")
    p.add_argument("--gff", default=env_default(
        "TFBS_GFF", "input_data/00_Reference_Genome/paramormyrops_ncbi_0.4.renamed.gff"),
        help="Annotation GFF; feature type CDS is excluded (default: $TFBS_GFF)")
    p.add_argument("--out_dir", default=str(ROOT / "output_data/15_TFBS_Scan"), help="Output directory")
    p.add_argument("--n_singletons", type=int, default=300000, help="Pool k = 1 units (default 300000)")
    p.add_argument("--n_chains", type=int, default=20000, help="Pool chains per k >= 2 (default 20000)")
    p.add_argument("--kmax_cap", type=int, default=4, help="KMAX = min(this, largest peak k) (default 4)")
    p.add_argument("--af_min", type=float, default=0.15, help="Keep pool SNPs with af > this (default 0.15)")
    p.add_argument("--peak_buffer", type=int, default=100000, help="Exclusion around peaks (default 100 kb)")
    p.add_argument("--seed", type=int, default=20261013, help="Random seed (default 20261013)")
    return p.parse_args()


def main():
    args = parse_args()
    t0 = time.time()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = IndexedFasta(args.ref)
    rng = np.random.default_rng(args.seed)

    # Peak units
    snps = pd.read_csv(args.polarized, sep="\t")
    peak_units, cds_dropped = build_peak_units(snps)
    fetch = fasta.fetch
    for u in peak_units:
        getter = lambda s, e, c=u.chr: fetch(c, s, e)  # noqa: E731
        getter.length = fasta.length(u.chr)
        annotate(u, getter)
    kmax = min(args.kmax_cap, max(u.k for u in peak_units))
    print("Peak units (cluster size k: n units)")
    for set_ in ("extended", "core"):
        for locus in LOCI:
            us = [u for u in peak_units if u.locus == locus and u.set == set_]
            dist = Counter(u.k for u in us)
            n_n = sum("N" in u.seq0 for u in us)
            print(f"  {set_:<8} {locus:<8} units={len(us):>3} snps={sum(u.k for u in us):>3}  "
                  + " ".join(f"k{k}:{dist[k]}" for k in sorted(dist))
                  + (f"  (units with reference N: {n_n})" if n_n else ""))
    print("CDS peak SNPs dropped: " + ", ".join(f"{l}={cds_dropped[l]}" for l in LOCI))
    print(f"KMAX = min({args.kmax_cap}, largest peak k = {max(u.k for u in peak_units)}) = {kmax}")

    # Pool
    t1 = time.time()
    df, counts = load_pool(args.assoc, args.af_min, fasta)
    peaks = read_peaks_bed(args.peaks_bed, args.peak_buffer)
    cds = read_cds_intervals(args.gff)
    pool, n_runs = filter_pool(df, peaks, cds, fasta, counts)
    del df
    print(f"\nPool filtering (read + filter {time.time() - t1:.0f} s)")
    for name, n in counts.items():
        print(f"  {name:<34}{n:>10,}")
    pool_units, stats = sample_pool_units(pool, n_runs, fasta, kmax, args.n_singletons, args.n_chains, rng)
    print("Pool units (k: chains of consecutive pool SNPs with gaps <= F / of those, span +-F in "
          "chromosome and N-free / sampled)")
    for k, n_chain, n_elig, take in stats:
        print(f"  k={k}: {n_chain:>10,} {n_elig:>10,} {take:>8,}")

    # Sequences and features for pool units, one chromosome at a time
    by_chr = {}
    for u in pool_units:
        by_chr.setdefault(u.chr, []).append(u)
    for chrom in fasta.order:
        if chrom not in by_chr:
            continue
        arr = fasta.chrom_array(chrom)
        getter = lambda s, e, a=arr: a[s - 1:e].tobytes().decode("ascii")  # noqa: E731
        getter.length = len(arr)
        for u in by_chr[chrom]:
            annotate(u, getter)
        del arr

    units = peak_units + pool_units
    for i, u in enumerate(units):
        u.idx = i
    tsv_path, fa_path = out_dir / "units.tsv", out_dir / "units.fa"
    write_outputs(units, tsv_path, fa_path)
    print(f"\nWrote {len(units):,} units ({len(peak_units)} peak, {len(pool_units):,} pool) to {tsv_path}")

    t2 = time.time()
    n_seq, n_bp = validate(tsv_path, fa_path, IndexedFasta(args.ref))
    print(f"Validation passed ({time.time() - t2:.0f} s): every unit and sequence re-read from disk "
          "and checked against the reference")
    print(f"units.fa: {n_seq:,} sequences, {n_bp:,} bp "
          f"(peak {2 * len(peak_units)} seqs; pool "
          + ", ".join(f"k{k}={sum(1 for u in pool_units if u.k == k) * 2 ** k:,}" for k in range(1, kmax + 1))
          + ")")
    print(f"Total runtime {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
