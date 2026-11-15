#!/usr/bin/env python3
"""24b_empirical_clustering.py -- empirical clustering test for Bambomo Peaks 1, 4 and 5.

Qualifying SNPs (from 24a): Bambomo - TP1 BP-derived delta >= the peak's index delta (Peak 4 0.8084, also used for
Peak 1; Peak 5 0.7836). For each qualifying SNP, k = other qualifying SNPs within +/- R. SNPs whose +/- R
neighbourhoods overlap are merged into blocks; block score = max k, focal SNP = the SNP with that k (lowest
position on ties). Peak score = the score of the block(s) overlapping the peak's extended interval. Null = blocks
clear of all six extended peaks +/- 250 kb.

Null sets: unmatched; primary match (SNP density within +/- 50 kb of the focal SNP within a factor 1.25, callable
fraction of its 1-Mb block within +/- 0.10; widened once to 1.5 / 0.15 if < 20 blocks); secondary = primary plus the
same tercile of TP1 LD (mean r^2 over 10-50 kb within +/- 100 kb of the focal SNP, ld23_lib, layout units).
p = (1 + null blocks scoring >= peak) / (1 + null blocks). R = 50 kb (primary) and 175 kb.

Usage (repo root): 24b_empirical_clustering.py <vcf> <24a_qualifying.tsv> <24a_density_1kb.tsv> <out_prefix>
"""
import bisect, csv, importlib.util, os, sys
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ld23_lib as L
spec = importlib.util.spec_from_file_location("obs23a", os.path.join(HERE, "23a_observed_ld_decay.py"))
obs23a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(obs23a)

BLOCKS = "output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree/sfs/A1w.block_scan.tsv"
PEAKS = "input_data/10_Introgression/peaks.bed"
TP1 = "output_data/11_Drift_vs_Selection/simulation/23_inputs/A1w_TP1.txt"
EPS = 1e-6                      # 24a writes delta to 6 decimals; distinct deltas differ by >= ~8e-5
THR = {"peak4": 124 / 134 - 11 / 94, "peak5": 105 / 134}
TESTS = [(1, "peak4"), (4, "peak4"), (5, "peak5")]
RADII = [50_000, 175_000]
PAD = 250_000
DENS_HALF = 50_000
LD_HALF = 100_000
LD_BINS = slice(5, 25)          # 10-50 kb
MATCH = [(1.25, 0.10), (1.5, 0.15)]
MIN_MATCHED = 20
MIN_READ = 19                   # below this, 1 / (n + 1) > 0.05 and the set is not read


def load(vcf, qual_fn, dens_fn):
    qual = defaultdict(list)
    for line in open(qual_fn):
        ch, p, d = line.split()
        qual[ch].append((int(p), float(d)))
    dens = defaultdict(dict)
    for line in open(dens_fn):
        ch, b, n = line.split()
        dens[ch][int(b)] = dens[ch].get(int(b), 0) + int(n)
    frac = {}
    with open(BLOCKS) as fh:
        next(fh)
        for line in fh:
            ch, b, ls = line.split("\t")[:3]
            frac[(ch, int(b))] = int(ls) / 1e6
    peaks = {}
    for line in open(PEAKS):
        f = line.split()
        peaks[int(f[3])] = (f"chr{f[0]}", int(f[1]), int(f[2]))
    return qual, dens, frac, peaks


def blocks_for(qual, thr, R):
    out = []
    for ch, rows in qual.items():
        pos = np.array(sorted(p for p, d in rows if d >= thr - EPS), dtype=np.int64)
        if pos.size == 0:
            continue
        k = np.searchsorted(pos, pos + R, "right") - np.searchsorted(pos, pos - R, "left") - 1
        start = 0
        for i in range(1, pos.size + 1):
            if i == pos.size or pos[i] - pos[i - 1] > 2 * R:
                seg = slice(start, i)
                j = start + int(np.argmax(k[seg]))
                out.append(dict(chrom=ch, start=int(pos[start]), end=int(pos[i - 1]), n_qual=i - start,
                                score=int(k[j]), focal=int(pos[j])))
                start = i
    return out


def overlaps(b, ch, s, e):
    return b["chrom"] == ch and b["start"] <= e and b["end"] >= s


def density(dens, ch, focal):
    lo, hi = (focal - DENS_HALF) // 1000, (focal + DENS_HALF) // 1000
    return sum(dens[ch].get(b, 0) for b in range(lo, hi + 1))


ld_cache = {}


def tp1_ld(vcf, ch, focal, f):
    key = (ch, focal)
    if key not in ld_cache:
        pos, G = obs23a.genotypes(vcf, TP1, ch, max(0, focal - LD_HALF), focal + LD_HALF)
        if pos.size < 2:
            ld_cache[key] = float("nan")
        else:
            M = (G >= 0).astype(np.int8)
            X = np.where(G >= 0, G, 0).astype(np.int8)
            keep = L.filter_snps(X, M)
            n, t = L.binned_r2(pos[keep], X[keep], M[keep], f)
            ld_cache[key] = t[LD_BINS].sum() / n[LD_BINS].sum() if n[LD_BINS].sum() else float("nan")
    return ld_cache[key]


def main():
    vcf, qual_fn, dens_fn, out = sys.argv[1:5]
    qual, dens, frac, peaks = load(vcf, qual_fn, dens_fn)
    block_rows, test_rows = [], []
    for R in RADII:
        for thr_name in ("peak4", "peak5"):
            bl = blocks_for(qual, THR[thr_name], R)
            for b in bl:
                b["f"] = frac.get((b["chrom"], b["focal"] // 1_000_000), float("nan"))
                b["density"] = density(dens, b["chrom"], b["focal"])
                b["peak"] = next((pk for pk, (ch, s, e) in peaks.items() if overlaps(b, ch, s, e)), "")
                b["near_peak"] = any(overlaps(b, ch, s - PAD, e + PAD) for ch, s, e in peaks.values())
            null = [b for b in bl if not b["near_peak"]]
            for b in null + [b for b in bl if b["peak"] != ""]:
                b["tp1_ld"] = tp1_ld(vcf, b["chrom"], b["focal"], b["f"])
            ld_null = np.array([b["tp1_ld"] for b in null])
            cuts = np.nanpercentile(ld_null, [100 / 3, 200 / 3])
            for b in bl:
                b["ld_tercile"] = int(np.searchsorted(cuts, b["tp1_ld"])) if "tp1_ld" in b and np.isfinite(b["tp1_ld"]) else ""
                block_rows.append(dict(radius_kb=R // 1000, threshold=thr_name, **{k: b.get(k, "") for k in (
                    "chrom", "start", "end", "n_qual", "score", "focal", "density", "f", "tp1_ld", "ld_tercile",
                    "peak", "near_peak")}))
            for pk, t in TESTS:
                if t != thr_name:
                    continue
                pb = [b for b in bl if b["peak"] == pk]
                top = max(pb, key=lambda b: b["score"])
                widened = 0
                for fac, df in MATCH:
                    prim = [b for b in null if top["density"] / fac <= b["density"] <= top["density"] * fac
                            and abs(b["f"] - top["f"]) <= df]
                    if len(prim) >= MIN_MATCHED:
                        break
                    widened = 1
                sec = [b for b in prim if b["ld_tercile"] == top["ld_tercile"]]
                for name, s in (("unmatched", null), ("matched", prim), ("matched_ld", sec)):
                    n_ge = sum(b["score"] >= top["score"] for b in s)
                    test_rows.append(dict(peak=pk, threshold=thr_name, delta_min=round(THR[thr_name], 4),
                                          radius_kb=R // 1000, null_set=name, peak_score=top["score"],
                                          peak_focal=f"{top['chrom']}:{top['focal']}", peak_density=top["density"],
                                          peak_f=round(float(top["f"]), 3),
                                          peak_tp1_ld=round(float(top["tp1_ld"]), 4),
                                          peak_ld_tercile=top["ld_tercile"], n_null=len(s), n_null_ge=n_ge,
                                          p=(1 + n_ge) / (1 + len(s)), p_floor=1 / (1 + len(s)),
                                          readable=int(len(s) >= MIN_READ),
                                          null_max=max((b["score"] for b in s), default=""),
                                          widened=widened if name != "unmatched" else ""))
                    print(test_rows[-1], flush=True)
    for fn, rows in ((out + "_blocks.tsv", block_rows), (out + "_test.tsv", test_rows)):
        with open(fn, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()
