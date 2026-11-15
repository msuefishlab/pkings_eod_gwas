#!/usr/bin/env python3
"""ext_aggregate.py -- point estimates and chromosome-bootstrap 95% bands of the panel, per chromosome set.

Sets: genome (chr1..25), calibration, heldout (chrom_split.tsv). Within a set, chromosomes are taken in sorted
(lexicographic) name order as in 23e; 1,000 bootstrap resamples of chromosomes with replacement, drawn as in 23e
(rng = numpy.random.default_rng(20261030), one rng.integers(0, n, n) per resample), a fresh rng for every
(sample-or-pair, set), so all statistics of one unit share the same resamples. Ratio statistics are recomputed
from resampled sums; counts are resampled sums. Interchromosomal r^2 uses the chromosome-pair matrices: for
resample weights w, value = w'Sw / w'Nw (zero diagonal: copies of one chromosome are never paired).
95% band = 2.5 / 97.5 percentiles. Same function serves simulated replicates (pass their segment tables).

Usage: ext_aggregate.py <seg_dir> <chrom_split.tsv> <samples_def.tsv> <pairs_def.tsv> <out_panel.tsv> <out_curves.tsv>
"""
import csv, os, sys
from collections import defaultdict
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel_lib as P

NBOOT, SEED, NBIN = 1000, 20261030, 50


def boot_idx(n):
    rng = np.random.default_rng(SEED)
    return np.vstack([rng.integers(0, n, n) for _ in range(NBOOT)])


def load_segments(seg_dir):
    """returns sums[(kind, name)][chrom] -> scalar or array(NBIN)"""
    sums = defaultdict(lambda: defaultdict(float))
    arrs = defaultdict(lambda: defaultdict(lambda: np.zeros(NBIN)))
    for fn in sorted(os.listdir(seg_dir)):
        if not (fn.startswith("seg.") and fn.endswith(".tsv")):
            continue
        for line in open(os.path.join(seg_dir, fn)):
            ch, b, kind, name, i, v = line.rstrip("\n").split("\t")
            if i == "":
                sums[(kind, name)][ch] += float(v)
            else:
                arrs[(kind, name)][ch][int(i)] += float(v)
    return sums, arrs


def load_inter(seg_dir, names):
    drawn = {n: {} for n in names}
    for fn in os.listdir(seg_dir):
        if fn.startswith("inter.") and fn.endswith(".npz") and ".tmp." not in fn:
            c = int(fn.split(".")[1][3:])
            z = np.load(os.path.join(seg_dir, fn))
            for n in names:
                drawn[n][c] = tuple(z[f"{n}__{k}"].astype(float) for k in ("X", "M", "Xp", "Mp"))
    return drawn


def band(x):
    return np.nanpercentile(x, 2.5, axis=0), np.nanpercentile(x, 97.5, axis=0)


def main():
    seg_dir, split_fn, sdef, pdef, out_panel, out_curves = sys.argv[1:7]
    split = {r["chrom"]: r["set"] for r in csv.DictReader(open(split_fn), delimiter="\t")}
    samples = [l.split("\t")[0] for l in open(sdef) if l.strip()]
    pairs = [tuple(l.rstrip("\n").split("\t")) for l in open(pdef) if l.strip()]
    sums, arrs = load_segments(seg_dir)
    allc = sorted(split)
    SETS = {"genome": allc, "calibration": sorted(c for c in allc if split[c] == "calibration"),
            "heldout": sorted(c for c in allc if split[c] == "heldout")}
    Lc = sums[("L", "")]
    print("callable bp per set:", {s: int(sum(Lc[c] for c in cs)) for s, cs in SETS.items()}, flush=True)
    drawn = load_inter(seg_dir, samples)
    inter = {}
    for n in samples:
        cs, S, N, Sp, Np = P.interchrom_matrices(drawn[n])
        inter[n] = ({f"chr{c}": i for i, c in enumerate(cs)}, S, N, Sp, Np)
    prow, crow = [], []

    def add(stat, unit, sname, val, bt, support):
        lo, hi = band(bt)
        prow.append(dict(statistic=stat, sample=unit, set=sname, value=val, lo=float(lo), hi=float(hi),
                         n_support=support))

    for sname, cs in SETS.items():
        idx = boot_idx(len(cs))
        Lv = np.array([Lc[c] for c in cs])
        Lb = Lv[idx].sum(1)
        for n in samples:
            idx = boot_idx(len(cs))
            Spi = np.array([sums[("S_pi", n)][c] for c in cs])
            add("pi", n, sname, Spi.sum() / Lv.sum(), Spi[idx].sum(1) / Lv[idx].sum(1), int(Lv.sum()))
            seg = np.array([sums[("n_segsites", n)][c] for c in cs])
            add("segsites_per_kb", n, sname, 1e3 * seg.sum() / Lv.sum(), 1e3 * seg[idx].sum(1) / Lv[idx].sum(1),
                int(seg.sum()))
            Nn = np.vstack([arrs[("ld_cnt", n)][c] for c in cs]); T = np.vstack([arrs[("ld_tot", n)][c] for c in cs])
            obs = T.sum(0) / Nn.sum(0)
            bt = np.vstack([T[i].sum(0) / Nn[i].sum(0) for i in idx])
            lo, hi = band(bt)
            for b in range(NBIN):
                crow.append(dict(kind="observed", pop=n, model="", r_over_mu="", bin=b, mid_kb=(b + 0.5) * 2,
                                 r2=obs[b], lo=lo[b], hi=hi[b], set=sname, n_pairs=int(Nn[:, b].sum())))
            far = T[:, P.FAR].sum(1); farn = Nn[:, P.FAR].sum(1)
            add("ld_r2_mean_20_100kb", n, sname, far.sum() / farn.sum(), far[idx].sum(1) / farn[idx].sum(1),
                int(farn.sum()))
            add("ld_r2_0_2kb", n, sname, obs[0], bt[:, 0], int(Nn[:, 0].sum()))
            pos, S, N, Sp, Np = inter[n]
            k = np.array([pos[c] for c in cs])
            S_, N_, Sp_, Np_ = (m[np.ix_(k, k)] for m in (S, N, Sp, Np))
            W = np.vstack([np.bincount(i, minlength=len(cs)) for i in idx]).astype(float)
            q = lambda Mx: np.einsum("bi,ij,bj->b", W, Mx, W)
            v_o, v_p = S_.sum() / N_.sum(), Sp_.sum() / Np_.sum()
            b_o, b_p = q(S_) / q(N_), q(Sp_) / q(Np_)
            npair = int(N_.sum() / 2)
            add("interchrom_r2", n, sname, v_o, b_o, npair)
            add("interchrom_r2_perm", n, sname, v_p, b_p, int(Np_.sum() / 2))
            add("interchrom_r2_excess", n, sname, v_o - v_p, b_o - b_p, npair)
        for bp, tp in pairs:
            key = f"{bp}|{tp}"
            idx = boot_idx(len(cs))
            g = lambda kind: np.array([sums[(kind, key)][c] for c in cs])
            dxy = g("S_dxy")
            add("dxy", key, sname, dxy.sum() / Lv.sum(), dxy[idx].sum(1) / Lv[idx].sum(1), int(Lv.sum()))
            fn_, fd_ = g("S_fst_num"), g("S_fst_den")
            add("fst_hudson", key, sname, fn_.sum() / fd_.sum(), fn_[idx].sum(1) / fd_[idx].sum(1), int(Lv.sum()))
            for tag in ("", "_24arule"):
                ns = g("delta_n_sites" + tag)
                add("delta_n_sites" + tag, key, sname, ns.sum(), ns[idx].sum(1), int(ns.sum()))
                for t in P.THRESH:
                    x = g(f"delta_ge_{t:g}" + tag)
                    add(f"delta_ge_{t:g}_count" + tag, key, sname, x.sum(), x[idx].sum(1), int(ns.sum()))
                    if tag == "":
                        add(f"delta_ge_{t:g}_per_100Mb", key, sname, 1e8 * x.sum() / Lv.sum(),
                            1e8 * x[idx].sum(1) / Lv[idx].sum(1), int(Lv.sum()))
    for fn, rows in ((out_panel, prow), (out_curves, crow)):
        with open(fn + ".tmp", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
            w.writeheader(); w.writerows(rows)
        os.replace(fn + ".tmp", fn)
    print("rows", len(prow), len(crow))


if __name__ == "__main__":
    main()
