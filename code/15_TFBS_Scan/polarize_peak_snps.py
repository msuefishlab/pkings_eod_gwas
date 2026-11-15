#!/usr/bin/env python3
"""
Stage A1 of the TFBS scan (specification: code/15_TFBS_Scan/README.md).

Polarize the extended GWAS peak SNPs (-log10 p > 3) that fall inside the three
TFBS windows: for each SNP, decide which allele is BP-associated and which is
TP-associated from ALT allele frequencies among BP and TP fish. No individual
haplotype is used downstream; the BP and TP haplotypes are the reference with
the BP or TP allele substituted at the peak SNPs.

Per SNP (one output row each):
  af_alt_bp, af_alt_tp  ALT frequency among called BP / TP fish (phased and
                        unphased GTs treated alike; missing GTs skipped)
  bp_allele             ALT if af_alt_bp > af_alt_tp, otherwise REF
  delta_af              |af_alt_bp - af_alt_tp| (= BP-allele frequency in BP
                        minus that in TP; always >= 0)
  af_ref                GEMMA `af` from the metadata, copied as README.md
                        specifies. NOTE: this `af` tracks the ALT (allele0)
                        frequency, not REF (allele1): corr = -0.98 with the REF
                        frequency among BP+TP fish, |af - ALT freq| <= 0.11.
                        Peak and pool carry the same quantity, so af-bin
                        matching is unaffected; do not read it as REF freq.
  in_cds                position inside a GFF feature of type CDS (only)
  core                  log_p > 5
  ref_ok                metadata REF equals the reference FASTA base

TFBS windows are 0-based half-open; a SNP is in a window if start < pos <= end.
Phenotypes: fish_data column 1 = sample id, column 7 = Phenotype; only `BP` and
`TP` fish are used.

Checks: VCF REF/ALT vs metadata allele1/allele0 and vs the FASTA (mismatches
are reported; an unresolvable mismatch is fatal), and the metadata `af`
against the REF and ALT frequencies among the BP+TP fish in the VCF (prints
both correlations and says which allele `af` tracks).

Usage:
  python3 code/15_TFBS_Scan/polarize_peak_snps.py \\
    [--metadata $TFBS_SNP_METADATA] [--vcf $TFBS_VCF] [--ref $TFBS_REF] \\
    [--gff $TFBS_GFF] [--fish_data input_data/01_Terra/data_model/fish_data_2023.txt] \\
    [--out output_data/15_TFBS_Scan/peak_snps_polarized.tsv]

Defaults come from the environment variables set by config.env when they
are defined, otherwise from the repo-relative paths that file sets.
"""

import argparse
import bisect
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

WINDOWS = {  # TFBS windows, 0-based half-open
    "sptbn4b": ("chr17", 11642886, 11742887),
    "cntnap5": ("chr16", 10839999, 11040000),
    "itih3": ("chr8", 15736955, 15836956),
}
LOCI = list(WINDOWS)
CORE_LOGP = 5.0

OUT_COLUMNS = [
    "locus", "chr", "pos", "snp_id", "ref", "alt", "log_p", "af_ref", "n_bp", "n_tp",
    "af_alt_bp", "af_alt_tp", "delta_af", "bp_allele", "tp_allele", "bp_is_ref",
    "in_cds", "core", "ref_ok",
]


# ── Shared helpers (also imported by build_units.py) ──────────────────────────
def env_default(var: str, rel_path: str) -> str:
    """Path from environment variable `var` if set, else ROOT/rel_path."""
    val = os.environ.get(var)
    return val if val else str(ROOT / rel_path)


class IndexedFasta:
    """Minimal random-access FASTA reader driven by a samtools .fai index."""

    def __init__(self, path: str):
        self.path = path
        self.index = {}
        self.order = []
        with open(path + ".fai") as fh:
            for line in fh:
                name, length, offset, linebases, linewidth = line.split("\t")[:5]
                self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))
                self.order.append(name)
        self.fh = open(path, "rb")

    def length(self, chrom: str) -> int:
        return self.index[chrom][0]

    def _file_offset(self, chrom: str, pos0: int) -> int:
        _, offset, lb, lw = self.index[chrom]
        return offset + (pos0 // lb) * lw + pos0 % lb

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Uppercase sequence of chrom:start-end (1-based, inclusive)."""
        length = self.index[chrom][0]
        if start < 1 or end > length or start > end:
            raise ValueError(f"bad interval {chrom}:{start}-{end} (length {length})")
        f0 = self._file_offset(chrom, start - 1)
        f1 = self._file_offset(chrom, end - 1) + 1
        self.fh.seek(f0)
        seq = self.fh.read(f1 - f0).replace(b"\n", b"").replace(b"\r", b"")
        if len(seq) != end - start + 1:
            raise RuntimeError(f"FASTA read length mismatch at {chrom}:{start}-{end}")
        return seq.decode("ascii").upper()

    def chrom_array(self, chrom: str) -> np.ndarray:
        """Whole chromosome as an uppercase uint8 array (index 0 = position 1)."""
        length = self.index[chrom][0]
        f0 = self._file_offset(chrom, 0)
        f1 = self._file_offset(chrom, length - 1) + 1
        self.fh.seek(f0)
        seq = self.fh.read(f1 - f0).replace(b"\n", b"").replace(b"\r", b"").upper()
        if len(seq) != length:
            raise RuntimeError(f"FASTA length mismatch for {chrom}: {len(seq)} vs {length}")
        return np.frombuffer(seq, dtype=np.uint8)


def read_cds_intervals(gff_path: str, chromosomes=None) -> dict:
    """
    {chrom: (starts, ends)} of merged CDS intervals, 1-based inclusive, as
    sorted int64 arrays. Only GFF feature type `CDS` is used (UTRs and other
    exon parts are not excluded). `chromosomes` restricts the parse if given.
    """
    raw = {}
    with open(gff_path) as fh:
        for line in fh:
            if "\tCDS\t" not in line or line.startswith("#"):
                continue
            parts = line.split("\t", 5)
            if parts[2] != "CDS":
                continue
            chrom = parts[0]
            if chromosomes is not None and chrom not in chromosomes:
                continue
            raw.setdefault(chrom, []).append((int(parts[3]), int(parts[4])))
    merged = {}
    for chrom, ivs in raw.items():
        ivs.sort()
        starts, ends = [], []
        for s, e in ivs:
            if starts and s <= ends[-1] + 1:
                ends[-1] = max(ends[-1], e)
            else:
                starts.append(s)
                ends.append(e)
        merged[chrom] = (np.asarray(starts, dtype=np.int64), np.asarray(ends, dtype=np.int64))
    return merged


def in_intervals(chrom: str, pos: int, intervals: dict) -> bool:
    """True if 1-based `pos` lies in one of the merged intervals of `chrom`."""
    if chrom not in intervals:
        return False
    starts, ends = intervals[chrom]
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and pos <= ends[i]


def in_intervals_vec(pos: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Vectorized in_intervals for one chromosome."""
    if len(starts) == 0:
        return np.zeros(len(pos), dtype=bool)
    i = np.searchsorted(starts, pos, side="right") - 1
    ok = i >= 0
    out = np.zeros(len(pos), dtype=bool)
    out[ok] = pos[ok] <= ends[i[ok]]
    return out


# ── Stage A1 ──────────────────────────────────────────────────────────────────
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Polarize peak SNPs in the TFBS windows (TFBS stage A1).")
    p.add_argument("--metadata", default=env_default(
        "TFBS_SNP_METADATA", "input_data/15_TFBS_Scan/PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS.txt"),
        help="Peak SNP table (default: $TFBS_SNP_METADATA)")
    p.add_argument("--vcf", default=env_default(
        "TFBS_VCF", "output_data/15_TFBS_Scan/PKINGS_ALL_WOB_EXCLUDED_SNPS_IN_PEAKS_TOP250.vcf.gz"),
        help="Indexed VCF with genotypes at the peak SNPs (default: $TFBS_VCF)")
    p.add_argument("--ref", default=env_default(
        "TFBS_REF", "input_data/00_Reference_Genome/"
        "jordan-mic4273-mb-hirise-20xx3__06-13-2024__final_assembly.fasta.renamed.fa"),
        help="Reference FASTA with .fai (default: $TFBS_REF)")
    p.add_argument("--gff", default=env_default(
        "TFBS_GFF", "input_data/00_Reference_Genome/paramormyrops_ncbi_0.4.renamed.gff"),
        help="Annotation GFF; feature type CDS is used (default: $TFBS_GFF)")
    p.add_argument("--fish_data", default=str(ROOT / "input_data/01_Terra/data_model/fish_data_2023.txt"),
                   help="Sample metadata; column 1 = sample id, column 7 = Phenotype")
    p.add_argument("--out", default=str(ROOT / "output_data/15_TFBS_Scan/peak_snps_polarized.tsv"),
                   help="Output TSV")
    p.add_argument("--bcftools", default=shutil.which("bcftools") or "/opt/homebrew/bin/bcftools",
                   help="bcftools executable")
    return p.parse_args()


def load_window_snps(metadata_path: str) -> pd.DataFrame:
    """Metadata rows inside the TFBS windows, with a `locus` column."""
    meta = pd.read_csv(metadata_path, sep="\t")
    parts = []
    for locus, (chrom, start, end) in WINDOWS.items():
        sub = meta[(meta["chr"] == chrom) & (meta["ps"] > start) & (meta["ps"] <= end)].copy()
        sub.insert(0, "locus", locus)
        parts.append(sub.sort_values("ps"))
    snps = pd.concat(parts, ignore_index=True)
    if snps.duplicated(["chr", "ps"]).any():
        sys.exit("ERROR: duplicate chr/ps in peak SNP metadata")
    return snps


def load_phenotypes(fish_data_path: str) -> dict:
    """{sample_id: 'BP' | 'TP'} from fish_data columns 1 and 7."""
    fish = pd.read_csv(fish_data_path, sep="\t", dtype=str)
    ids, pheno = fish.iloc[:, 0], fish.iloc[:, 6]
    return {s: p for s, p in zip(ids, pheno) if p in ("BP", "TP")}


def query_vcf(bcftools: str, vcf: str) -> tuple:
    """(sample list, {(chr, pos): [(ref, alt, [gt, ...]), ...]}) for the three windows."""
    samples = subprocess.run([bcftools, "query", "-l", vcf], check=True,
                             capture_output=True, text=True).stdout.split()
    regions = ",".join(f"{c}:{s + 1}-{e}" for c, s, e in WINDOWS.values())
    out = subprocess.run(
        [bcftools, "query", "-r", regions, "-f", r"%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n", vcf],
        check=True, capture_output=True, text=True).stdout
    records = {}
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) != 4 + len(samples):
            sys.exit(f"ERROR: unexpected bcftools field count at {f[0]}:{f[1]}")
        records.setdefault((f[0], int(f[1])), []).append((f[2], f[3], f[4:]))
    return samples, records


def alt_counts(gts: list, idx: np.ndarray) -> tuple:
    """(ALT allele count, called allele count, called individuals) over samples `idx`."""
    n_alt = n_called_alleles = n_ind = 0
    for i in idx:
        alleles = gts[i].replace("|", "/").split("/")
        if "." in alleles:
            continue
        if any(a not in ("0", "1") for a in alleles):
            raise ValueError(f"non-biallelic genotype {gts[i]}")
        n_ind += 1
        n_called_alleles += len(alleles)
        n_alt += sum(a == "1" for a in alleles)
    return n_alt, n_called_alleles, n_ind


def polarize(snps: pd.DataFrame, samples: list, records: dict, pheno: dict,
             fasta: IndexedFasta, cds: dict) -> tuple:
    """Build the output table; returns (table, list of problem messages)."""
    bp_idx = np.array([i for i, s in enumerate(samples) if pheno.get(s) == "BP"])
    tp_idx = np.array([i for i, s in enumerate(samples) if pheno.get(s) == "TP"])
    print(f"VCF samples: {len(samples)}; BP fish in VCF: {len(bp_idx)}; TP fish in VCF: {len(tp_idx)}")

    rows, problems, vcf_ref_freq = [], [], []
    for r in snps.itertuples(index=False):
        chrom, pos, ref, alt = r.chr, int(r.ps), str(r.allele1), str(r.allele0)
        recs = records.get((chrom, pos), [])
        if len(recs) != 1:
            problems.append(f"{r.rs}: {len(recs)} VCF records at {chrom}:{pos}")
            if not recs:
                continue
        vref, valt, gts = recs[0]
        flip = False
        if (vref, valt) != (ref, alt):
            if (vref, valt) == (alt, ref):
                flip = True
                problems.append(f"{r.rs}: VCF REF/ALT {vref}/{valt} swapped vs metadata {ref}/{alt}")
            else:
                problems.append(f"{r.rs}: VCF REF/ALT {vref}/{valt} != metadata {ref}/{alt}; SNP skipped")
                continue
        a_bp, n_bp_al, n_bp = alt_counts(gts, bp_idx)
        a_tp, n_tp_al, n_tp = alt_counts(gts, tp_idx)
        if n_bp_al == 0 or n_tp_al == 0:
            problems.append(f"{r.rs}: no called BP or TP genotypes; SNP skipped")
            continue
        af_alt_bp, af_alt_tp = a_bp / n_bp_al, a_tp / n_tp_al
        if flip:  # VCF ALT is metadata REF
            af_alt_bp, af_alt_tp = 1.0 - af_alt_bp, 1.0 - af_alt_tp
        vcf_ref_freq.append(1.0 - (a_bp + a_tp) / (n_bp_al + n_tp_al) if not flip
                            else (a_bp + a_tp) / (n_bp_al + n_tp_al))
        bp_is_ref = not (af_alt_bp > af_alt_tp)
        fasta_base = fasta.fetch(chrom, pos, pos)
        rows.append({
            "locus": r.locus, "chr": chrom, "pos": pos, "snp_id": r.rs, "ref": ref, "alt": alt,
            "log_p": float(r.log_p), "af_ref": float(r.af), "n_bp": n_bp, "n_tp": n_tp,
            "af_alt_bp": round(af_alt_bp, 6), "af_alt_tp": round(af_alt_tp, 6),
            "delta_af": round(abs(af_alt_bp - af_alt_tp), 6),
            "bp_allele": ref if bp_is_ref else alt, "tp_allele": alt if bp_is_ref else ref,
            "bp_is_ref": bp_is_ref, "in_cds": in_intervals(chrom, pos, cds),
            "core": float(r.log_p) > CORE_LOGP, "ref_ok": fasta_base == ref,
        })
        if fasta_base != ref:
            problems.append(f"{r.rs}: metadata REF {ref} != FASTA base {fasta_base}")
    table = pd.DataFrame(rows, columns=OUT_COLUMNS)
    if len(table) > 2:
        af, ref_freq = table["af_ref"].to_numpy(), np.asarray(vcf_ref_freq)
        rho = np.corrcoef(af, ref_freq)[0, 1]
        d_ref, d_alt = np.max(np.abs(af - ref_freq)), np.max(np.abs(af - (1.0 - ref_freq)))
        print(f"Check: corr(metadata af, VCF REF freq among BP+TP) = {rho:+.3f}; "
              f"max |af - REF freq| = {d_ref:.3f}, max |af - ALT freq| = {d_alt:.3f}")
        print("       -> metadata `af` tracks the " + ("REF (allele1)" if rho > 0 else "ALT (allele0)")
              + " frequency; written to af_ref unchanged, as README.md specifies")
    return table, problems


def print_summary(table: pd.DataFrame) -> None:
    print("\nlocus     n  n_in_cds  n_delta<0.1  frac_bp_is_ref  n_core  ref_ok_failures")
    for locus in LOCI:
        t = table[table["locus"] == locus]
        print(f"{locus:<8}{len(t):>4}{int(t['in_cds'].sum()):>10}{int((t['delta_af'] < 0.1).sum()):>13}"
              f"{t['bp_is_ref'].mean():>11.3f} ({int(t['bp_is_ref'].sum())}/{len(t)})"
              f"{int(t['core'].sum()):>6}{int((~t['ref_ok']).sum()):>10}")


def main():
    args = parse_args()
    snps = load_window_snps(args.metadata)
    print(f"Peak SNPs in TFBS windows: {len(snps)} "
          + ", ".join(f"{l}={int((snps['locus'] == l).sum())}" for l in LOCI))
    pheno = load_phenotypes(args.fish_data)
    samples, records = query_vcf(args.bcftools, args.vcf)
    fasta = IndexedFasta(args.ref)
    cds = read_cds_intervals(args.gff, {c for c, _, _ in WINDOWS.values()})
    table, problems = polarize(snps, samples, records, pheno, fasta, cds)

    print(f"\nProblems ({len(problems)}):" if problems else "\nProblems: none (VCF REF/ALT = metadata "
          "allele1/allele0 = FASTA at every SNP)")
    for msg in problems:
        print("  " + msg)
    print_summary(table)

    if len(table) != len(snps):
        sys.exit(f"ERROR: {len(snps) - len(table)} window SNPs could not be polarized; nothing written")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote {len(table)} rows to {args.out}")


if __name__ == "__main__":
    main()
