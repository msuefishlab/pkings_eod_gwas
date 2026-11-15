#!/usr/bin/env python3
"""ext_mask_split.py -- model extension, Stage 1: summary-statistic mask and calibration/held-out split.

Mask: exactly the callable sequence simulation/23a drew its segments from, i.e. every 1-Mb block
returned by 23a_observed_ld_decay.eligible_blocks() (imported unchanged): blocks of the drift test's A1w block scan
(fixed_differences_agree/sfs/A1w.block_scan.tsv) on chr1..chr25 with L_scan > 0, lying wholly inside their
chromosome, and not overlapping the chr4 inversion-like polymorphism or any of the six extended peaks
(input_data/10_Introgression/peaks.bed, cols 2-3) +/- 250 kb. Callable length of a block = L_scan (A1w callable
sites: AN >= projection in BAM and TP1, outgroup-polarized, ExcHet >= 0.001, outside A1w call deserts).
23a scored a random 200 of these blocks; the panel uses all of them.

Split: chromosomes sorted by masked callable length (descending; ties by chromosome number), paired
consecutively; per pair one rng.permutation(2) from numpy.random.default_rng(20261108), in sorted order;
perm[0] -> calibration, perm[1] -> heldout; the unpaired last chromosome -> calibration.

Usage: ext_mask_split.py <code_dir_with_23a> <block_scan.tsv> <peaks.bed> <chrom_lengths.tsv> <out_dir>
Writes mask_blocks.tsv, chrom_split.tsv, ld23a_draw.tsv (the 200 blocks 23a's seed selects).
"""
import importlib.util, os, sys
import numpy as np

code, block_scan, peaks, lens_fn, out = sys.argv[1:6]
spec = importlib.util.spec_from_file_location("a23", os.path.join(code, "23a_observed_ld_decay.py"))
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
lengths = {l.split()[0]: int(l.split()[1]) for l in open(lens_fn)}
el = A.eligible_blocks(block_scan, peaks, lengths)
os.makedirs(out, exist_ok=True)
with open(os.path.join(out, "mask_blocks.tsv"), "w") as fo:
    fo.write("chrom\tblock\tstart\tend\tf\tL_callable\n")
    for ch, b, f in el:
        fo.write(f"{ch}\t{b}\t{b * 1_000_000}\t{(b + 1) * 1_000_000}\t{f:.6f}\t{int(round(f * 1_000_000))}\n")
rng = np.random.default_rng(A.SEED)
draw = [el[i] for i in sorted(rng.choice(len(el), A.N_SEG, replace=False))]
with open(os.path.join(out, "ld23a_draw.tsv"), "w") as fo:
    fo.write("chrom\tblock\tf\n")
    for ch, b, f in draw:
        fo.write(f"{ch}\t{b}\t{f:.6f}\n")
L = {}
for ch, b, f in el:
    L[ch] = L.get(ch, 0) + int(round(f * 1_000_000))
chroms = sorted(L, key=lambda c: (-L[c], int(c[3:])))
assert len(chroms) == 25, chroms
srng = np.random.default_rng(20261108)
rows = []
for k in range(0, 24, 2):
    p = srng.permutation(2)
    pair = (chroms[k], chroms[k + 1])
    rows.append((pair[p[0]], "calibration", k // 2 + 1))
    rows.append((pair[p[1]], "heldout", k // 2 + 1))
rows.append((chroms[24], "calibration", 13))
nb = {}
for ch, b, f in el:
    nb[ch] = nb.get(ch, 0) + 1
with open(os.path.join(out, "chrom_split.tsv"), "w") as fo:
    fo.write("chrom\tset\tpair\trank\tL_callable_masked\tn_blocks\n")
    for ch, s, pr in rows:
        fo.write(f"{ch}\t{s}\t{pr}\t{chroms.index(ch) + 1}\t{L[ch]}\t{nb[ch]}\n")
for s in ("calibration", "heldout"):
    cs = [ch for ch, t, _ in rows if t == s]
    print(s, len(cs), "chromosomes", sum(L[c] for c in cs), "callable bp", sum(nb[c] for c in cs), "blocks")
print("genome", sum(L.values()), "callable bp", len(el), "blocks")
