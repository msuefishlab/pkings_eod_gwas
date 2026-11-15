#!/usr/bin/env python3
"""10_shared_segments.py -- Bambomo-like segments shared by the 2019 Apassa wobbles (post hoc).

Carrier at a marker: either haplotype with posterior P(BAM) >= 0.5 (hap_primary/post_all.npz, all markers).
Qualifying marker: all three 2019-018/021 wobbles carry, >= 8 of the 10 2019 wobbles carry, and <= 2 of the 9
2009 Apassa TP fish outside the panels carry. Qualifying markers within 500 kb merge into loci.
Null: 1,000 independent circular rotations of each fish's carrier track along the marker index (keeps coverage
and tract structure, breaks alignment between fish). p = (1 + rotations with >= observed loci) / 1001.
Also reports, for each locus, its span, marker count and carrier counts; and the local-ancestry posterior of
every Apassa fish at chr7:13,388,648 (nearest marker), the existing Apassa wobble lead.

Writes <out_dir>/results/shared_segments.tsv and results/shared_segments_summary.txt.

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/10_shared_segments.py <out_dir>
"""
import sys
import numpy as np

SEED = 20261103
NROT = 1000


def loci_from(q, chrom, pos, gap=500_000):
    idx = np.flatnonzero(q)
    out = []
    for i in idx:
        if out and chrom[i] == chrom[out[-1][-1]] and pos[i] - pos[out[-1][-1]] <= gap:
            out[-1].append(i)
        else:
            out.append([i])
    return out


def main():
    O = sys.argv[1]
    info = {}
    for line in open(f"{O}/samples.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "sample":
            info[f[0]] = (f[2], f[4], set(f[5].split(",")))
    z = np.load(f"{O}/hap_primary/post_all.npz")
    post, chrom, pos, hs = z["post"], z["chrom"], z["pos"], z["hap_sample"]
    fish = list(dict.fromkeys(hs))
    carr = np.stack([post[:, hs == s].max(1) >= 0.5 for s in fish], 1)
    W = [i for i, s in enumerate(fish) if info[s][0] == "APA2019_Wobble"]
    S18 = [i for i in W if info[fish[i]][1] != "2019-017"]
    C = [i for i, s in enumerate(fish) if info[s][0] == "APA2009_TP"]
    assert len(W) == 10 and len(S18) == 3 and len(C) == 9

    def qualify(X):
        return X[:, S18].all(1) & (X[:, W].sum(1) >= 8) & (X[:, C].sum(1) <= 2)

    q = qualify(carr)
    loci = loci_from(q, chrom, pos)
    rng = np.random.default_rng(SEED)
    n = carr.shape[0]
    cols = W + C
    null = []
    for _ in range(NROT):
        X = carr.copy()
        for j in cols:
            X[:, j] = np.roll(carr[:, j], rng.integers(0, n))
        null.append(len(loci_from(qualify(X), chrom, pos)))
    null = np.array(null)
    p = (1 + (null >= len(loci)).sum()) / (1 + NROT)
    rows = ["locus\tchrom\tstart\tend\tn_markers\tW19_carriers_max\tS18_carriers\tC09_carriers_min\tpeak"]
    for k, L in enumerate(loci, 1):
        L = np.array(L)
        rows.append(f"{k}\tchr{chrom[L[0]]}\t{pos[L].min()}\t{pos[L].max()}\t{L.size}\t{carr[L][:, W].sum(1).max()}\t"
                    f"{carr[L][:, S18].sum(1).min()}\t{carr[L][:, C].sum(1).min()}\t{int(z['peak'][L].max())}")
    open(f"{O}/results/shared_segments.tsv", "w").write("\n".join(rows) + "\n")
    lines = [f"qualifying markers {q.sum()} of {n}; loci {len(loci)}; rotation null median {np.median(null):.0f}, "
             f"95th pct {np.percentile(null, 95):.0f}, max {null.max()}; p = {p:.4f}"]
    # per-marker carrier frequencies, for context
    lines.append(f"genome-wide carrier fraction per marker: 2019 wobbles {carr[:, W].mean():.3f}, "
                 f"2019-018/021 {carr[:, S18].mean():.3f}, 2009 TP {carr[:, C].mean():.3f}")
    allthree = carr[:, S18].all(1)
    lines.append(f"markers where all three 2019-018/021 carry: {allthree.sum()} ({allthree.mean():.3f})")
    # chr7 lead
    s7 = np.flatnonzero(chrom == 7)
    j = s7[np.argmin(np.abs(pos[s7] - 13_388_648))]
    lines.append(f"chr7:13,388,648 nearest marker chr7:{pos[j]} ({pos[j] - 13_388_648:+d} bp); P(BAM) h0/h1 per fish:")
    for s in fish:
        if info[s][0].startswith("APA"):
            v = post[j, hs == s]
            lines.append(f"  {s}\t{info[s][0]}\t{info[s][1]}\t{v[0]:.3f}/{v[1]:.3f}")
    open(f"{O}/results/shared_segments_summary.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(rows)); print("\n".join(lines))


if __name__ == "__main__":
    main()
