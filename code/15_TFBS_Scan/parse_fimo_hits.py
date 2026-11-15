#!/usr/bin/env python3
"""
Stream FIMO --text output and keep the best hit per (unit, mask, motif) that covers a unit SNP.

FIMO must run with --text --no-pgc, so start/stop are 1-based positions within each unit
sequence (README.md, stage F). Every kept hit is checked against the unit sequence
(matched_sequence must equal the substring, reverse-complemented on the minus strand);
any mismatch aborts, because it means the coordinates are wrong.

Two modes:
  parse:  fimo ... chunk.fa | parse_fimo_hits.py parse --units units.tsv --fasta chunk.fa \\
              --motifs JASPAR.meme --out chunk.npz
  merge:  parse_fimo_hits.py merge --motifs JASPAR.meme --out hits.npz chunk_*.npz
"""

import argparse
import sys

import numpy as np
import pandas as pd

COMP = str.maketrans("ACGTN", "TGCAN")


def parse_args():
    p = argparse.ArgumentParser(description="Parse FIMO --text hits covering unit SNPs.")
    sub = p.add_subparsers(dest="mode", required=True)
    a = sub.add_parser("parse")
    a.add_argument("--units", required=True)
    a.add_argument("--fasta", required=True)
    a.add_argument("--motifs", required=True)
    a.add_argument("--out", required=True)
    m = sub.add_parser("merge")
    m.add_argument("--motifs", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("chunks", nargs="+")
    return p.parse_args()


def read_motif_order(path):
    ids, names = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("MOTIF"):
                parts = line.split()
                ids.append(parts[1])
                names.append(parts[2] if len(parts) > 2 else parts[1])
    return ids, names


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip().upper())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def parse(args):
    motif_ids, _ = read_motif_order(args.motifs)
    motif_index = {m: i for i, m in enumerate(motif_ids)}
    seqs = read_fasta(args.fasta)
    units = pd.read_csv(args.units, sep="\t", usecols=["unit_idx", "snp_off"], dtype={"snp_off": str})
    offsets = dict(zip(units["unit_idx"].to_numpy(),
                       (tuple(int(x) for x in s.split(",")) for s in units["snp_off"])))

    best = {}
    n_lines = n_cover = 0
    for line in sys.stdin:
        if not line or line[0] == "#" or line.startswith("motif_id"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        n_lines += 1
        name = f[2]
        start, stop = int(f[3]), int(f[4])
        unit = int(name[1:8])
        if not any(start <= o <= stop for o in offsets[unit]):
            continue
        n_cover += 1
        frag = seqs[name][start - 1:stop]
        if f[5] == "-":
            frag = frag.translate(COMP)[::-1]
        if frag != f[9].upper():
            sys.exit(f"[ERROR] matched_sequence mismatch: {line.strip()} vs {frag}")
        key = (unit, int(name[10:]), motif_index[f[0]])
        mlp = -np.log10(float(f[7]))
        if mlp > best.get(key, -1.0):
            best[key] = mlp

    keys = np.array(list(best.keys()), dtype=np.int64).reshape(-1, 3)
    np.savez_compressed(
        args.out,
        unit_idx=keys[:, 0].astype(np.int32),
        mask=keys[:, 1].astype(np.int16),
        motif_idx=keys[:, 2].astype(np.int16),
        mlog10p=np.array(list(best.values()), dtype=np.float32),
    )
    print(f"[INFO] {args.fasta}: {n_lines} hits, {n_cover} cover a SNP, "
          f"{len(best)} (unit, mask, motif) kept", file=sys.stderr)


def merge(args):
    motif_ids, motif_names = read_motif_order(args.motifs)
    parts = [np.load(c) for c in args.chunks]
    out = {k: np.concatenate([p[k] for p in parts]) for k in ["unit_idx", "mask", "motif_idx", "mlog10p"]}
    key = (out["unit_idx"].astype(np.int64) << 24) | (out["mask"].astype(np.int64) << 12) | out["motif_idx"]
    if len(np.unique(key)) != len(key):
        sys.exit("[ERROR] duplicate (unit, mask, motif) keys across chunks")
    np.savez_compressed(args.out, **out,
                        motif_ids=np.array(motif_ids), motif_names=np.array(motif_names))
    print(f"[INFO] merged {len(args.chunks)} chunks: {len(key)} rows, "
          f"{len(np.unique(out['unit_idx']))} units with hits -> {args.out}")


def main():
    args = parse_args()
    parse(args) if args.mode == "parse" else merge(args)


if __name__ == "__main__":
    main()
