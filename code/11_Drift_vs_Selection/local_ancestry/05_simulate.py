#!/usr/bin/env python3
"""05_simulate.py -- simulated admixed fish from held-out haplotypes, and the calibration they give.

Mosaics are spliced from held-out haplotypes only (4 clean Apassa TP fish = 8 haplotypes, 10 Bambomo 2009
fish = 20 haplotypes; none in any panel of the primary marker set) at the primary panel's genome-wide markers.
Along each chromosome, at r = 1.4 cM/Mb, segments are exponential with rate T per Morgan; at each breakpoint the
ancestry is redrawn (BAM with probability m), and the donor is redrawn only when the ancestry changes
so a same-ancestry stretch copies one real haplotype. The chain starts at stationarity.
No-tract control (T = inf): every marker drawn independently. Grid:
T in {1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000, 3000, inf} x m in {0.05, 0.10, 0.15, 0.20, 0.25},
plus m = 0 (T label 0); 10 diploid fish per cell. Seed 20261103.

Each fish gets the same statistics as the real fish (04): phased HMM fit (m_hat, T_hat), q_m projection score,
called tracts and F_long (fraction of analysed length in tracts >= 2 Mb); and the diploid HMM fit on the summed
genotypes (m_hat_dip, T_hat_dip).

Calibration: leave-one-out matching of every simulated fish against the cell statistics (la_lib.match)
gives the matched T; T is recovered in a cell if >= 8 of 10 fish have a best-match T within [T/3, 3T]; the power
range is the contiguous set of recovered T at m = 0.15.

Writes <out_dir>/sim/: fish.tsv, recovery.tsv (per cell and model), power.txt.

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/05_simulate.py <out_dir> [workers]
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import la_lib as L

T_SIM = [1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000, 3000, np.inf]
M_SIM = [0.05, 0.10, 0.15, 0.20, 0.25]
NREP = 10
SEED = 20261103
R_SIM = 1.4e-8


def read_roles(out):
    roles = {}
    for line in open(f"{out}/samples.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "sample":
            roles[f[0]] = set(f[5].split(","))
    return roles


def simulate_hap(rng, T, m, pos, seg, HA, HB):
    """Returns alleles (n,), true ancestry at markers (n,), true BAM intervals clipped to segments."""
    n = pos.size
    allele = np.empty(n, np.int8); anc = np.empty(n, np.int8)
    intervals = []
    for s in np.unique(seg):
        ii = np.flatnonzero(seg == s)
        p = pos[ii]
        if np.isinf(T):
            a = rng.random(ii.size) < m
            dA = rng.integers(0, HA.shape[1], ii.size); dB = rng.integers(0, HB.shape[1], ii.size)
            anc[ii] = a
            allele[ii] = np.where(a, HB[ii, dB], HA[ii, dA])
            continue
        x0, x1 = float(p[0]), float(p[-1])
        if T == 0:                                   # m = 0 cell: one Apassa donor per segment
            bps = np.array([x0, x1 + 1])
        else:
            bps = [x0]
            while bps[-1] <= x1:
                bps.append(bps[-1] + rng.exponential(1.0 / T) / R_SIM)
            bps = np.array(bps)
        k = bps.size
        a = rng.random(k) < m
        dB, dA = np.empty(k, int), np.empty(k, int)
        for j in range(k):                           # donor kept while the ancestry stays the same
            if j and a[j] == a[j - 1]:
                dB[j], dA[j] = dB[j - 1], dA[j - 1]
            else:
                dB[j], dA[j] = rng.integers(0, HB.shape[1]), rng.integers(0, HA.shape[1])
        which = np.searchsorted(bps, p, side="right") - 1
        anc[ii] = a[which]
        allele[ii] = np.where(a[which], HB[ii, dB[which]], HA[ii, dA[which]])
        cur = None
        for j in range(k):
            lo, hi = bps[j], min(bps[j + 1] if j + 1 < k else x1, x1)
            if hi <= lo:
                continue
            if a[j]:
                cur = [lo, hi] if cur is None else [cur[0], hi]
            elif cur is not None:
                intervals.append(tuple(cur)); cur = None
        if cur is not None:
            intervals.append(tuple(cur))
    return allele, anc, intervals


def contiguous_best(lst):
    best, cur = [], []
    for T, ok in lst:
        cur = cur + [T] if ok else []
        if len(cur) > len(best):
            best = cur
    return best


def main():
    out = sys.argv[1]
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    od = f"{out}/sim"; os.makedirs(od, exist_ok=True)
    roles = read_roles(out)
    M = L.load_markers(out, "primary")
    gw = M["gw"].astype(bool)
    hs = M["hap_sample"]
    colA = np.flatnonzero([("holdout_APA" in roles.get(s, ())) for s in hs])
    colB = np.flatnonzero([("holdout_BAM" in roles.get(s, ())) for s in hs])
    assert colA.size == 8 and colB.size == 20
    chrom, pos = M["chrom"][gw], M["pos"][gw]
    pBe, pTe = M["pBe"][gw], M["pTe"][gw]
    HA, HB = M["H"][gw][:, colA], M["H"][gw][:, colB]
    seg = L.segments(chrom, pos, L.excluded_regions())
    alen = L.analysed_length(pos, seg)

    rng = np.random.default_rng(SEED)
    cells = [(0, 0.0)] + [(T, m) for m in M_SIM for T in T_SIM]
    hap_obs, truth, meta = [], [], []
    for T, m in cells:
        for rep in range(NREP):
            for _ in range(2):
                al, an, iv = simulate_hap(rng, T, m, pos, seg, HA, HB)
                hap_obs.append(al); truth.append((an, iv))
            meta.append((T, m, rep))
    obs = np.stack(hap_obs, 1)
    geno = obs[:, 0::2] + obs[:, 1::2]
    print(f"simulated {len(meta)} fish on {pos.size} markers", flush=True)

    ones = np.ones(pos.size, bool)
    Sh = L.grid_surfaces("hap", L.per_chrom_split(chrom, pos, obs, pBe, pTe, ones), workers)
    Sh = Sh[..., 0::2] + Sh[..., 1::2]
    fits_h = L.fit_from_surfaces(Sh, B=200); del Sh
    Sd = L.grid_surfaces("dip", L.per_chrom_split(chrom, pos, geno, pBe, pTe, ones), workers)
    fits_d = L.fit_from_surfaces(Sd, B=200); del Sd
    print("fitted", flush=True)

    lam = np.repeat([f["T"] * L.R_FIT for f in fits_h], 2)
    mm = np.repeat([f["m"] for f in fits_h], 2)
    post = np.empty(obs.shape)
    for c in range(1, 26):
        s = chrom == c
        post[s] = L.hap_posterior(obs[s], pos[s], pBe[s], pTe[s], lam, mm)

    rows = []
    for i, (T, m, rep) in enumerate(meta):
        fh, fd = fits_h[i], fits_d[i]
        sm, _ = L.summarise_fish(post[:, 2 * i:2 * i + 2], pos, chrom, seg, fh["T"], fh["m"])
        ivs = truth[2 * i][1] + truth[2 * i + 1][1]
        tl = np.array([b - a for a, b in ivs]) if ivs else np.zeros(0)
        an = np.stack([truth[2 * i][0], truth[2 * i + 1][0]], 1).astype(bool)
        call = post[:, 2 * i:2 * i + 2] >= 0.5
        rows.append(dict(T=T, m=m, rep=rep, true_frac_markers=float(an.mean()),
                         true_frac_len=float(tl.sum() / (2 * alen)), true_n=int(tl.size),
                         true_mean_tract=float(tl.mean()) if tl.size else np.nan,
                         true_F_long=float(tl[tl >= L.LONG_BP].sum() / (2 * alen)) if tl.size else 0.0,
                         q_m=L.projection_q(obs[:, 2 * i:2 * i + 2], pBe, pTe),
                         m_hat=fh["m"], T_hat=fh["T"], logT=float(np.log10(fh["T"])),
                         frac_post=sm["frac_post"], frac_tract=sm["frac_tract"], F_long=sm["F_long"],
                         n_tracts=sm["n_tracts"], median_tract=sm["med_tract"], longest=sm["max_tract"],
                         sens=float((call & an).sum() / max(an.sum(), 1)),
                         prec=float((call & an).sum() / max(call.sum(), 1)),
                         m_hat_dip=fd["m"], T_hat_dip=fd["T"], logT_dip=float(np.log10(fd["T"]))))

    # leave-one-out matching
    for model, keys in L.MATCH_KEYS.items():
        cs = L.cell_stats(rows, keys)
        for i, r in enumerate(rows):
            c = (r["T"], r["m"])
            same = [j for j, x in enumerate(rows) if (x["T"], x["m"]) == c]
            best, compat, _ = L.match(r, cs, keys, exclude=(c, same.index(i)))
            r[f"match_T_{model}"] = best[0]; r[f"match_m_{model}"] = best[1]
            Ts = [x[0] for x in compat]
            r[f"match_T_lo_{model}"] = min(Ts) if Ts else np.nan
            r[f"match_T_hi_{model}"] = max(Ts) if Ts else np.nan
    keys = list(rows[0])
    with open(f"{od}/fish.tsv", "w") as fh_:
        fh_.write("\t".join(keys) + "\n")
        for r in rows:
            fh_.write("\t".join(f"{r[k]:.6g}" if isinstance(r[k], (float, np.floating)) else str(r[k]) for k in keys) + "\n")

    rec = ["model\tT\tm\tn\tmedian_q_m\tmedian_m_hat\tmedian_T_hat\tmedian_F_long\tmedian_true_F_long\t"
           "median_match_T\tmatch_T_p10\tmatch_T_p90\tn_match_within_3x\trecovered\tmedian_match_m\t"
           "n_interval_covers_T\tmedian_sens\tmedian_prec"]
    power = {}
    for model in L.MATCH_KEYS:
        mk = "m_hat" if model == "phased" else "m_hat_dip"
        tk = "T_hat" if model == "phased" else "T_hat_dip"
        for T, m in cells:
            rr = [r for r in rows if r["T"] == T and r["m"] == m]
            mt = np.array([r[f"match_T_{model}"] for r in rr], float)
            lo = np.array([r[f"match_T_lo_{model}"] for r in rr], float)
            hi = np.array([r[f"match_T_hi_{model}"] for r in rr], float)
            if T in (0,) or np.isinf(T):
                within = int(np.sum(mt == T)); ok = within >= 8
            else:
                within = int(np.sum((mt >= T / 3) & (mt <= 3 * T))); ok = within >= 8
            cover = int(np.sum((lo <= T) & (hi >= T)))
            if m == 0.15:
                power.setdefault(model, []).append((T, ok))
            rec.append(f"{model}\t{T}\t{m}\t{len(rr)}\t{np.median([r['q_m'] for r in rr]):.4f}\t"
                       f"{np.median([r[mk] for r in rr]):.3f}\t{np.median([r[tk] for r in rr]):.1f}\t"
                       f"{np.median([r['F_long'] for r in rr]):.4f}\t{np.median([r['true_F_long'] for r in rr]):.4f}\t"
                       f"{np.median(mt):.4g}\t{np.percentile(mt, 10):.4g}\t{np.percentile(mt, 90):.4g}\t{within}\t"
                       f"{int(ok)}\t{np.median([r[f'match_m_{model}'] for r in rr]):.3f}\t{cover}\t"
                       f"{np.median([r['sens'] for r in rr]):.3f}\t{np.median([r['prec'] for r in rr]):.3f}")
    open(f"{od}/recovery.tsv", "w").write("\n".join(rec) + "\n")
    lines = []
    for model, lst in power.items():
        okT = [T for T, ok in lst if ok]
        best = [T for T in contiguous_best(lst) if np.isfinite(T)]
        if best:
            lo, hi = best[0], best[-1]
            lines.append(f"{model}\tm=0.15\trecovered_T={okT}\tpower_range_r1.4={lo}-{hi}\t"
                         f"r1.0={lo * 1.4:.3g}-{hi * 1.4:.3g}\tr3.0={lo * 1.4 / 3:.3g}-{hi * 1.4 / 3:.3g}")
        else:
            lines.append(f"{model}\tm=0.15\trecovered_T={okT}\tpower_range=none")
    open(f"{od}/power.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(rec)); print("\n".join(lines))


if __name__ == "__main__":
    main()
