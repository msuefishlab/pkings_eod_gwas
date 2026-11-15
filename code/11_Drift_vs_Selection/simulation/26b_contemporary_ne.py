#!/usr/bin/env python3
"""26b_contemporary_ne.py -- contemporary Ne in Bambomo and Apassa.

Reads the per-chromosome genotype tables of 26a_extract_temporal.sh. Its Ne estimates are not reported in
the paper; 26c_ancestry_check.py imports its loader (load) and site filters (MAF, CALL_T).

Temporal method, per population (2009 vs 2019): Nei-Tajima Fc as a ratio of sums over loci,
Fc = sum (x - y)^2 / sum ((x + y)/2 - x y), x and y the 2009 and 2019 alt-allele frequencies.
Waples (1989) sampling plan II: drift = Fc - 1/(2 S0) - 1/(2 St), Ne = t / (2 drift), with S the mean
number of called individuals and t = 10 years / generation time (1 and 2 years reported).
95% CI from a delete-one-chromosome jackknife of the drift term. Sites: call rate >= 0.8 in both
year-samples and pooled MAF >= 0.05.
  Control 1: 200 random splits of the 2009 sample into two groups of the 2019 sample's size and the
  rest; the split's Fc minus its analytic sampling term should centre on 0.
  Control 2: Spearman correlation, across sites kept in both populations, of the 2009->2019 frequency
  change in Bambomo and in Apassa (about 0 under independent drift; positive suggests a shared
  batch artefact or gene flow).

LD method, per year-sample: up to 300 evenly spaced SNPs per chromosome (MAF >= 0.05, call rate >= 0.9
in that sample); r^2 = genotype correlation over jointly called individuals (>= 10), pairs on
different chromosomes only. Sampling term from 3 permutations of the individuals of the second
chromosome of each pair; r2' = r2_obs - r2_perm. Ne by Waples (2006): S >= 30,
Ne = (1/3 + sqrt(1/9 - 2.76 r2'))/(2 r2'); S < 30, Ne = (0.308 + sqrt(0.308^2 - 2.08 r2'))/(2 r2').
95% CI by bootstrap over chromosomes (1,000 replicates), on r2' and then converted to Ne.
Each sample is run twice: call rate >= 0.9 and complete data (at call rate 0.9, r2_obs < r2_perm in every
bootstrap replicate, a bias from missingness that is correlated within fish).

Usage (repo root): 26b_contemporary_ne.py <vcf> <26_temporal dir> <out_prefix>
"""
import gzip, os, subprocess, sys
import numpy as np
from scipy.stats import spearmanr

YEARS = 10
GEN = [1, 2]
MAF, CALL_T, CALL_LD, MIN_JOINT = 0.05, 0.8, 0.9, 10
N_LD, N_PERM, N_BOOT, N_SPLIT = 300, 3, 1000, 200
rng = np.random.default_rng(20261103)


def load(vcf, d):
    order = subprocess.run(f"bcftools view -h -S {d}/samples.txt {vcf} 2>/dev/null | grep '^#CHROM'",
                           shell=True, capture_output=True, text=True, check=True).stdout.split()[9:]
    grp = dict(line.split() for line in open(f"{d}/groups.tsv"))
    labels = np.array([grp[s] for s in order])
    code = {"0/0": 0, "0/1": 1, "1/0": 1, "1/1": 2, "0|0": 0, "0|1": 1, "1|0": 1, "1|1": 2}
    data = {}
    for c in range(1, 26):
        ch = f"chr{c}"
        pos, rows = [], []
        with gzip.open(f"{d}/geno/{ch}.tsv.gz", "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                gts = f[2:]
                if not any(g in ("0/1", "1/1", "0|1", "1|0", "1|1", "1/0") for g in gts):
                    continue
                pos.append(int(f[1]))
                rows.append([code.get(g, -1) for g in gts])
        data[ch] = (np.array(pos), np.array(rows, dtype=np.int8))
        print(f"  {ch}: {len(pos)} polymorphic records", flush=True)
    return labels, data


def counts(G, cols):
    sub = G[:, cols]
    ok = sub >= 0
    return np.where(ok, sub, 0).sum(1).astype(float), ok.sum(1).astype(float)


def fc_parts(a0, n0, a1, n1):
    x, y = a0 / (2 * n0), a1 / (2 * n1)
    return ((x - y) ** 2).sum(), ((x + y) / 2 - x * y).sum(), n0.sum(), n1.sum(), x.size


def drift_from(parts):
    num, den, s0, s1, k = (sum(p[i] for p in parts) for i in range(5))
    S0, S1 = s0 / k, s1 / k
    fc = num / den
    return fc, S0, S1, fc - 1 / (2 * S0) - 1 / (2 * S1)


def temporal(labels, data, pop):
    c0, c1 = np.flatnonzero(labels == f"{pop}2009"), np.flatnonzero(labels == f"{pop}2019")
    per_chr, keep_all, freq = {}, {}, {}
    for ch, (pos, G) in data.items():
        a0, n0 = counts(G, c0)
        a1, n1 = counts(G, c1)
        with np.errstate(invalid="ignore", divide="ignore"):
            pool = (a0 + a1) / (2 * (n0 + n1))
            keep = (n0 >= CALL_T * c0.size) & (n1 >= CALL_T * c1.size) & (pool >= MAF) & (pool <= 1 - MAF)
        per_chr[ch] = fc_parts(a0[keep], n0[keep], a1[keep], n1[keep])
        keep_all[ch] = keep
        freq[ch] = dict(zip(pos[keep], (a1[keep] / (2 * n1[keep]) - a0[keep] / (2 * n0[keep]))))
    fc, S0, S1, dr = drift_from(per_chr.values())
    jk = np.array([drift_from([v for k, v in per_chr.items() if k != ch])[3] for ch in per_chr])
    k = len(jk)
    se = np.sqrt((k - 1) / k * ((jk - jk.mean()) ** 2).sum())
    n_snps = sum(v[4] for v in per_chr.values())
    rows = []
    for g in GEN:
        t = YEARS / g
        ne = t / (2 * dr) if dr > 0 else float("inf")
        lo_d, hi_d = dr - 1.96 * se, dr + 1.96 * se
        rows.append(dict(pop=pop, generation_years=g, t_generations=t, n_snps=n_snps, S0=round(S0, 2),
                         St=round(S1, 2), Fc=fc, sampling_term=1 / (2 * S0) + 1 / (2 * S1), drift=dr,
                         drift_se=se, Ne=ne, Ne_lo=t / (2 * hi_d) if hi_d > 0 else float("inf"),
                         Ne_hi=t / (2 * lo_d) if lo_d > 0 else float("inf")))
    # control 1: pseudo-temporal splits of the 2009 sample
    split = []
    for _ in range(N_SPLIT):
        perm = rng.permutation(c0)
        ga, gb = perm[c1.size:], perm[:c1.size]
        parts = []
        for ch, (pos, G) in data.items():
            kp = keep_all[ch]
            a0, n0 = counts(G[kp], ga)
            a1, n1 = counts(G[kp], gb)
            ok = (n0 > 0) & (n1 > 0)
            parts.append(fc_parts(a0[ok], n0[ok], a1[ok], n1[ok]))
        split.append(drift_from(parts)[3])
    return rows, np.array(split), freq


def ld_cross(XA, MA, XB, MB):
    n = MA @ MB.T
    sx, sy = XA @ MB.T, MA @ XB.T
    sxx, syy, sxy = (XA ** 2) @ MB.T, MA @ (XB ** 2).T, XA @ XB.T
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = sxy / n - sx * sy / n ** 2
        vx, vy = sxx / n - (sx / n) ** 2, syy / n - (sy / n) ** 2
        r2 = cov ** 2 / (vx * vy)
    ok = (n >= MIN_JOINT) & np.isfinite(r2) & (vx > 0) & (vy > 0)
    return r2[ok].sum(), ok.sum(), n[ok].sum()


def waples(r2p, S):
    if r2p <= 0:
        return float("inf")
    if S >= 30:
        disc = 1 / 9 - 2.76 * r2p
        return (1 / 3 + np.sqrt(disc)) / (2 * r2p) if disc >= 0 else float("nan")
    disc = 0.308 ** 2 - 2.08 * r2p
    return (0.308 + np.sqrt(disc)) / (2 * r2p) if disc >= 0 else float("nan")


def ld_ne(labels, data, group, call_ld=CALL_LD):
    cols = np.flatnonzero(labels == group)
    blocks = []
    for ch, (pos, G) in data.items():
        sub = G[:, cols]
        ok = sub >= 0
        n = ok.sum(1)
        with np.errstate(invalid="ignore", divide="ignore"):
            p = np.where(ok, sub, 0).sum(1) / (2 * n)
        keep = np.flatnonzero((n >= call_ld * cols.size) & (p >= MAF) & (p <= 1 - MAF))
        if keep.size > N_LD:
            keep = keep[np.linspace(0, keep.size - 1, N_LD).astype(int)]
        X = np.where(ok[keep], sub[keep], 0).astype(float)
        blocks.append((X, ok[keep].astype(float)))
    C = len(blocks)
    So, Sp, Np, Nj = (np.zeros((C, C)) for _ in range(4))
    for a in range(C):
        for b in range(a + 1, C):
            XA, MA = blocks[a]; XB, MB = blocks[b]
            s, k, nj = ld_cross(XA, MA, XB, MB)
            So[a, b], Np[a, b], Nj[a, b] = s, k, nj
            ps = 0.0
            for _ in range(N_PERM):
                pm = rng.permutation(cols.size)
                s_p, k_p, _ = ld_cross(XA, MA, XB[:, pm], MB[:, pm])
                ps += s_p / k_p
            Sp[a, b] = ps / N_PERM * k
    def est(idx):
        so = sp = npairs = nj = 0.0
        for i in range(len(idx)):
            for j in range(i + 1, len(idx)):
                a, b = sorted((idx[i], idx[j]))
                if a == b:
                    continue
                so += So[a, b]; sp += Sp[a, b]; npairs += Np[a, b]; nj += Nj[a, b]
        r2o, r2p = so / npairs, sp / npairs
        S = nj / npairs
        return r2o, r2p, S, waples(r2o - r2p, S)
    r2o, r2p, S, ne = est(list(range(C)))
    # CI on r2' (finite even when Ne is unbounded), then converted: the upper r2' bound gives the lower Ne bound
    boots = []
    for _ in range(N_BOOT):
        bo, bp, _, _ = est(list(rng.integers(0, C, C)))
        boots.append(bo - bp)
    r2_lo, r2_hi = np.percentile(boots, [2.5, 97.5])
    return dict(group=group, call_rate_min=call_ld, n_individuals=cols.size, mean_joint_S=round(S, 2),
                snps_per_chrom=int(np.mean([b[0].shape[0] for b in blocks])), r2_obs=r2o, r2_perm=r2p,
                r2_drift=r2o - r2p, r2_drift_lo=r2_lo, r2_drift_hi=r2_hi, Ne=ne,
                Ne_lo=waples(r2_hi, S), Ne_hi=waples(r2_lo, S), small_sample=int(cols.size < 25))


def write(fn, rows):
    with open(fn, "w") as fh:
        fh.write("\t".join(rows[0].keys()) + "\n")
        for r in rows:
            fh.write("\t".join(str(v) for v in r.values()) + "\n")


def main():
    vcf, d, out = sys.argv[1:4]
    labels, data = load(vcf, d)
    print({g: int((labels == g).sum()) for g in sorted(set(labels))}, flush=True)
    trows, split_rows, freqs = [], [], {}
    for pop in ("BAM", "APA"):
        rows, split, freqs[pop] = temporal(labels, data, pop)
        trows += rows
        split_rows.append(dict(pop=pop, n_splits=split.size, mean_split_drift=split.mean(),
                               sd_split_drift=split.std(), q025=np.percentile(split, 2.5),
                               q975=np.percentile(split, 97.5), observed_drift=rows[0]["drift"],
                               share_splits_ge_observed=float(np.mean(split >= rows[0]["drift"]))))
        for r in rows:
            print(r, flush=True)
        print(split_rows[-1], flush=True)
    xs, ys = [], []
    for ch in freqs["BAM"]:
        common = set(freqs["BAM"][ch]) & set(freqs["APA"][ch])
        xs += [freqs["BAM"][ch][p] for p in common]
        ys += [freqs["APA"][ch][p] for p in common]
    rho, pv = spearmanr(xs, ys)
    batch = dict(n_common_sites=len(xs), spearman_rho=rho, p=pv,
                 pearson_r=float(np.corrcoef(xs, ys)[0, 1]))
    print("batch check", batch, flush=True)
    ## call_ld 1.0: at 0.9, r2_obs < r2_perm in every bootstrap replicate, because missingness
    ## is correlated within fish, which the permutation breaks, biasing the sampling term upward
    lrows = [ld_ne(labels, data, g, c) for c in (CALL_LD, 1.0)
             for g in ("BAM2009", "APA2009", "BAM2019", "APA2019")]
    for r in lrows:
        print(r, flush=True)
    write(out + "_temporal.tsv", trows)
    write(out + "_split_control.tsv", split_rows)
    write(out + "_batch_check.tsv", [batch])
    write(out + "_ld.tsv", lrows)


if __name__ == "__main__":
    main()
