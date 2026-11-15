#!/usr/bin/env python3
"""
Check that the short-flank unit scans equal a full-window scan at the three peaks.

For each locus, build the full-window BP and TP haplotypes (reference plus BP or TP alleles
at every extended-set peak SNP), padded by F = 32 bp on each side, and run
FIMO --text --no-pgc on them. Hits covering a peak SNP are mapped to that SNP's unit and
reduced to the best p per (unit, mask, motif). The result must equal the peak rows of
hits.npz exactly; any difference means the flank construction or the coordinate handling is wrong.

Usage:
  python3 code/15_TFBS_Scan/verify_flank_equivalence.py \\
    --units output_data/15_TFBS_Scan/units.tsv --hits output_data/15_TFBS_Scan/hits.npz \\
    --ref $TFBS_REF --motifs $JASPAR_MOTIFS --bfile output_data/15_TFBS_Scan/genome.bg
"""

import argparse
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

F = 32
FIMO = "fimo"  # MEME Suite 5.5.8, on the PATH
COMP = str.maketrans("ACGTN", "TGCAN")


def parse_args():
    p = argparse.ArgumentParser(description="Flank-scan vs full-window FIMO equivalence check.")
    p.add_argument("--units", required=True)
    p.add_argument("--hits", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--motifs", required=True)
    p.add_argument("--bfile", required=True)
    return p.parse_args()


def fetch(ref, fai, chrom, start1, end1):
    """Reference bases start1..end1 (1-based inclusive), uppercase."""
    length, offset, line_bases, line_width = fai[chrom]
    out = []
    with open(ref, "rb") as fh:
        p0, p1 = start1 - 1, end1
        first_line, last_line = p0 // line_bases, (p1 - 1) // line_bases
        fh.seek(offset + first_line * line_width)
        raw = fh.read((last_line - first_line + 1) * line_width).replace(b"\n", b"").replace(b"\r", b"")
        out = raw[p0 - first_line * line_bases: p1 - first_line * line_bases]
    return out.decode().upper()


def main():
    args = parse_args()
    fai = {}
    with open(args.ref + ".fai") as fh:
        for line in fh:
            f = line.split("\t")
            fai[f[0]] = tuple(int(x) for x in f[1:5])
    motif_ids = [l.split()[1] for l in open(args.motifs) if l.startswith("MOTIF")]
    midx = {m: i for i, m in enumerate(motif_ids)}

    units = pd.read_csv(args.units, sep="\t", dtype=str)
    peak = units[(units["source"] == "peak") & (units["set"] == "extended")]
    h = np.load(args.hits)
    sel = np.isin(h["unit_idx"], peak["unit_idx"].astype(int).to_numpy())
    flank = {(int(u), int(m), int(k)): round(float(p), 3)
             for u, m, k, p in zip(h["unit_idx"][sel], h["mask"][sel], h["motif_idx"][sel], h["mlog10p"][sel])}

    full = {}
    for locus, g in peak.groupby("locus"):
        chrom = g["chr"].iloc[0]
        snps = []  # (pos, unit, bit, ref, alt, bp_mask, k)
        for _, r in g.iterrows():
            pos = [int(x) for x in r["pos"].split(",")]
            refs, alts = r["ref"].split(","), r["alt"].split(",")
            for j, (p, a, b) in enumerate(zip(pos, refs, alts)):
                snps.append((p, int(r["unit_idx"]), j, a, b, int(r["bp_mask"]), int(r["k"])))
        lo, hi = min(s[0] for s in snps) - F, max(s[0] for s in snps) + F
        base = fetch(args.ref, fai, chrom, lo, hi)
        by_pos = {s[0]: s for s in snps}
        haps = {}
        for hap in ("BP", "TP"):
            seq = list(base)
            for p, u, j, a, b, bp_mask, k in snps:
                assert seq[p - lo] == a, f"{locus} {p}: ref mismatch"
                bp_alt = (bp_mask >> j) & 1
                use_alt = bp_alt if hap == "BP" else 1 - bp_alt
                seq[p - lo] = b if use_alt else a
            haps[hap] = "".join(seq)
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tf:
            for hap, s in haps.items():
                tf.write(f">{hap}_{locus}\n{s}\n")
            fa = tf.name
        res = subprocess.run([FIMO, "--text", "--no-pgc", "--thresh", "1e-3", "--bfile", args.bfile,
                              args.motifs, fa], capture_output=True, text=True, check=True)
        for line in res.stdout.splitlines():
            if line.startswith(("motif_id", "#")) or not line.strip():
                continue
            f = line.split("\t")
            hap = f[2].split("_")[0]
            s1, e1 = int(f[3]) + lo - 1, int(f[4]) + lo - 1  # genomic 1-based
            frag = haps[hap][int(f[3]) - 1:int(f[4])]
            if f[5] == "-":
                frag = frag.translate(COMP)[::-1]
            assert frag == f[9].upper(), f"matched_sequence mismatch in full-window scan: {line}"
            covered_units = {by_pos[p][1] for p in range(s1, e1 + 1) if p in by_pos}
            for u in covered_units:
                r = peak[peak["unit_idx"] == str(u)].iloc[0]
                bp_mask, k = int(r["bp_mask"]), int(r["k"])
                mask = bp_mask if hap == "BP" else bp_mask ^ ((1 << k) - 1)
                key = (u, mask, midx[f[0]])
                p = round(-np.log10(float(f[7])), 3)
                if p > full.get(key, -1):
                    full[key] = p

    only_flank = set(flank) - set(full)
    only_full = set(full) - set(flank)
    diff_p = [k for k in set(flank) & set(full) if abs(flank[k] - full[k]) > 1e-3]
    print(f"[CHECK] peak (unit, mask, motif) rows: flank={len(flank)} full={len(full)} "
          f"only_flank={len(only_flank)} only_full={len(only_full)} p_diff={len(diff_p)}")
    if only_flank or only_full or diff_p:
        for k in list(only_flank)[:5] + list(only_full)[:5] + diff_p[:5]:
            print("  ", k, flank.get(k), full.get(k))
        sys.exit("[FAIL] flank scan differs from full-window scan")
    print("[PASS] flank scan == full-window scan at all peak units")


if __name__ == "__main__":
    main()
