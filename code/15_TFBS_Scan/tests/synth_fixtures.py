"""Synthetic inputs in README.md format (units.tsv, hits.npz, families table, peak_snps_polarized.tsv) for
test_tfbs_null.py and for runtime checks. No real data are read.

Hit model: each unit gets Poisson(hit_rate) motif hits with a base -log10 p of 3 + Exp(ln 10) (10% above 4); every
haplotype (mask) of the unit gets base + noise * N(0, 1), kept when >= 3. Haplotypes of one unit are exchangeable, and
peak units follow the same model as pool units, so the matched null is exact unless a signal is planted.
"""
from pathlib import Path

import numpy as np
import pandas as pd

F = 32
BASES = np.array(list("ACGT"))
LOCUS_CHR = {"sptbn4b": "chr17", "cntnap5": "chr16", "itih3": "chr8"}
LOCUS_START = {"sptbn4b": 11642886, "cntnap5": 10839999, "itih3": 15736955}
NAMED_IDS = ["MA0463.3", "MA0731.1", "MA0102.5", "MA0466.4", "MA0836.3", "MA0837.3", "MA0838.1", "MA1636.2",
             "MA0019.2"]
NAMED_NAMES = ["BCL6", "BCL6B", "CEBPA", "CEBPB", "CEBPD", "CEBPE", "CEBPG", "CEBPG", "Ddit3::Cebpa"]


def motif_table(n_motifs=60, n_fam=10):
    ids = NAMED_IDS + [f"MS{i:04d}.1" for i in range(n_motifs - len(NAMED_IDS))]
    names = NAMED_NAMES + [f"TF{i}" for i in range(n_motifs - len(NAMED_IDS))]
    fams, classes = [], []
    for i in range(n_motifs):
        if i < 2:
            fams.append("More than 3 adjacent zinc fingers")
            classes.append("C2H2 zinc finger factors")
        elif i < 9:
            fams.append("C/EBP")
            classes.append("Basic leucine zipper factors (bZIP)")
        else:
            j = (i - 9) % n_fam
            fams.append(f"Fam{j:02d}")
            classes.append(f"Class{j % 3}")
    fams[-1] = "Fam00;Fam01"      # several families: the first is primary
    fams[-2] = ""                 # empty family: Unknown
    return pd.DataFrame(dict(motif_id=ids, motif_name=names, tf_class=classes, tf_family=fams))


def _gap_bin(pos):
    m = int(np.min(np.diff(pos)))
    return "1-4" if m <= 4 else "5-12" if m <= 12 else "13-32"


def _rand_alleles(rng, k):
    ref = rng.integers(0, 4, k)
    alt = (ref + rng.integers(1, 4, k)) % 4
    return BASES[ref], BASES[alt]


def build_unit_records(rng, n_pool1, n_chain, kmax, peak_k, n_core, chain_gap_max=F):
    """Return a DataFrame with one row per unit (lists in object columns), peak units first."""
    recs = []
    for locus, ks in peak_k.items():
        chrom, pos = LOCUS_CHR[locus], LOCUS_START[locus] + 100
        ext = []
        for i, k in enumerate(ks):
            pos += int(rng.integers(40, 400))
            p = [pos]
            for _ in range(k - 1):
                p.append(p[-1] + int(rng.integers(1, F + 1)))
            pos = p[-1]
            ref, alt = _rand_alleles(rng, k)
            core = i < n_core.get(locus, 0) and k == 1
            lp = list(rng.uniform(5.2, 8.0, k) if core else rng.uniform(3.1, 4.9, k))
            ext.append(dict(source="peak", locus=locus, set="extended", chr=chrom, k=k, pos=p, ref=list(ref),
                            alt=list(alt), gc=float(rng.uniform(0.3, 0.6)), cpg=list(rng.integers(0, 2, k)),
                            af=list(rng.uniform(0.15, 0.85, k)), log_p=lp, bp_mask=int(rng.integers(0, 2 ** k)),
                            core=core))
        recs += ext
        recs += [dict(r, set="core") for r in ext if r["core"]]
    # pool singletons
    pos = np.cumsum(rng.integers(200, 2000, n_pool1)) + 1000
    ref, alt = _rand_alleles(rng, n_pool1)
    gc = rng.uniform(0.3, 0.6, n_pool1)
    cpg = (rng.random(n_pool1) < 0.25).astype(int)
    af = rng.uniform(0.15, 0.85, n_pool1)
    for i in range(n_pool1):
        recs.append(dict(source="pool", locus="pool", set="pool", chr="chrP1", k=1, pos=[int(pos[i])], ref=[ref[i]],
                         alt=[alt[i]], gc=float(gc[i]), cpg=[int(cpg[i])], af=[float(af[i])], log_p=None,
                         bp_mask=None))
    # pool chains
    start = 1000
    for k in range(2, kmax + 1):
        for _ in range(n_chain):
            start += int(rng.integers(300, 3000))
            p = [start]
            for _ in range(k - 1):
                p.append(p[-1] + int(rng.integers(1, chain_gap_max + 1)))
            start = p[-1]
            r, a = _rand_alleles(rng, k)
            recs.append(dict(source="pool", locus="pool", set="pool", chr="chrP2", k=k, pos=p, ref=list(r),
                             alt=list(a), gc=float(rng.uniform(0.3, 0.6)), cpg=list(rng.integers(0, 2, k)),
                             af=list(rng.uniform(0.15, 0.85, k)), log_p=None, bp_mask=None))
    df = pd.DataFrame(recs)
    df.insert(0, "unit_idx", np.arange(len(df)))
    return df


def write_units(df, path):
    """Write unit records (list columns) as units.tsv in the README.md format."""
    j = lambda xs, fmt="{}": ",".join(fmt.format(x) for x in xs)  # noqa: E731
    out = pd.DataFrame(dict(
        unit_idx=df["unit_idx"], source=df["source"], locus=df["locus"], set=df["set"], chr=df["chr"], k=df["k"],
        pos=[j(p) for p in df["pos"]], ref=[j(x) for x in df["ref"]], alt=[j(x) for x in df["alt"]],
        snp_off=[j([F + 1 + (q - p[0]) for q in p]) for p in df["pos"]],
        seq_start=[p[0] - F for p in df["pos"]], seq_len=[(p[-1] - p[0]) + 2 * F + 1 for p in df["pos"]],
        gap_bin=[_gap_bin(p) if len(p) > 1 else "NA" for p in df["pos"]], gc=[f"{x:.4f}" for x in df["gc"]],
        cpg=[j(x) for x in df["cpg"]], af_ref=[j(x, "{:.3f}") for x in df["af"]],
        snp_id=[j([f"{c}_{q}" for q in p]) for c, p in zip(df["chr"], df["pos"])],
        log_p=[j(x, "{:.4f}") if isinstance(x, list) else "NA" for x in df["log_p"]],
        bp_mask=[str(int(x)) if s == "peak" else "NA" for s, x in zip(df["source"], df["bp_mask"])]))
    out.to_csv(path, sep="\t", index=False)


def write_polarized(df, path):
    rows = []
    for r in df[df["source"] == "peak"].itertuples():
        for jj in range(r.k):
            bit = (int(r.bp_mask) >> jj) & 1
            rows.append(dict(locus=r.locus, chr=r.chr, pos=r.pos[jj], snp_id=f"{r.chr}_{r.pos[jj]}", ref=r.ref[jj],
                             alt=r.alt[jj], log_p=r.log_p[jj], af_ref=r.af[jj], n_bp=100, n_tp=100,
                             af_alt_bp=0.6 if bit else 0.4, af_alt_tp=0.4 if bit else 0.6, delta_af=0.2,
                             bp_allele=r.alt[jj] if bit else r.ref[jj], tp_allele=r.ref[jj] if bit else r.alt[jj],
                             bp_is_ref=bit == 0, in_cds=False, core=r.log_p[jj] > 5, ref_ok=True))
    pd.DataFrame(rows).drop_duplicates(["chr", "pos"]).to_csv(path, sep="\t", index=False)


def gen_hits(rng, df, n_mot, hit_rate, noise):
    k = df["k"].to_numpy()
    peak = (df["source"] == "peak").to_numpy()
    bpm = np.array([int(x) if p else 0 for x, p in zip(df["bp_mask"], peak)], dtype=np.int64)
    n = len(df)
    nb = rng.poisson(hit_rate * k, n)   # hits scale with the number of SNPs in the unit
    bu = np.repeat(np.arange(n), nb)
    bm = rng.integers(0, n_mot, len(bu))
    key = np.unique(bu.astype(np.int64) * n_mot + bm)
    bu, bm = key // n_mot, key % n_mot
    bs = 3 + rng.exponential(1 / np.log(10), len(bu))
    nm = np.where(peak, 2, 2 ** k)
    rep = np.repeat(np.arange(len(bu)), nm[bu])
    u = bu[rep]
    starts = np.r_[0, np.cumsum(nm[bu])[:-1]]
    slot = np.arange(len(rep)) - np.repeat(starts, nm[bu])
    full = 2 ** k[u] - 1
    mask = np.where(peak[u], np.where(slot == 0, bpm[u], bpm[u] ^ full), slot)
    s = bs[rep] + noise * rng.standard_normal(len(rep))
    keep = s >= 3
    return dict(unit_idx=u[keep], mask=mask[keep], motif_idx=bm[rep][keep], mlog10p=s[keep])


def plant(h, df, motifs, locus, n, plant_motif, set_name="extended"):
    """Make n peak units at `locus` BP-only for `plant_motif`: clear `motifs` (the family / group) on the TP haplotype
    and on the BP haplotype, then add plant_motif at -log10 p = 6 on the BP haplotype."""
    ids = list(motifs["motif_id"])
    clear = np.array([ids.index(x) for x in plant_motif["clear"]])
    pm = ids.index(plant_motif["motif"])
    sel = df[(df["source"] == "peak") & (df["locus"] == locus) & (df["set"] == set_name)].head(n)
    if len(sel) < n:
        raise ValueError("not enough peak units to plant")
    targets = sel["unit_idx"].to_numpy()
    drop = np.isin(h["unit_idx"], targets) & np.isin(h["motif_idx"], clear)
    out = {key: v[~drop] for key, v in h.items()}
    bpm = np.array([int(x) for x in sel["bp_mask"]])
    out["unit_idx"] = np.r_[out["unit_idx"], targets]
    out["mask"] = np.r_[out["mask"], bpm]
    out["motif_idx"] = np.r_[out["motif_idx"], np.full(n, pm)]
    out["mlog10p"] = np.r_[out["mlog10p"], np.full(n, 6.0)]
    return out


def write_hits(h, motifs, path):
    np.savez(path, unit_idx=np.asarray(h["unit_idx"], np.int32), mask=np.asarray(h["mask"], np.int16),
             motif_idx=np.asarray(h["motif_idx"], np.int16), mlog10p=np.asarray(h["mlog10p"], np.float32),
             motif_ids=np.array(motifs["motif_id"], dtype=str), motif_names=np.array(motifs["motif_name"], dtype=str))


DEFAULT_PEAK_K = {"sptbn4b": [1] * 50 + [2, 2, 3],
                  "cntnap5": [1] * 60 + [2, 2, 3, 6],
                  "itih3": [1] * 45 + [2, 3, 4, 5]}
DEFAULT_CORE = {"sptbn4b": 25, "cntnap5": 25, "itih3": 25}


def make_fixture(d, seed=1, n_pool1=6000, n_chain=600, kmax=4, peak_k=None, n_core=None, n_motifs=60, n_fam=10,
                 hit_rate=10.0, noise=0.6, plant_spec=None, write_pol=True):
    """Write units.tsv, hits.npz, families.tsv (and peak_snps_polarized.tsv) into directory d; return their paths
    plus the unit records and the motif table."""
    d = Path(d)
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    peak_k = DEFAULT_PEAK_K if peak_k is None else peak_k
    n_core = DEFAULT_CORE if n_core is None else n_core
    df = build_unit_records(rng, n_pool1, n_chain, kmax, peak_k, n_core)
    motifs = motif_table(n_motifs, n_fam)
    h = gen_hits(rng, df, n_motifs, hit_rate, noise)
    if plant_spec:
        h = plant(h, df, motifs, **plant_spec)
    paths = dict(units=d / "units.tsv", hits=d / "hits.npz", families=d / "families.tsv",
                 polarized=d / "peak_snps_polarized.tsv" if write_pol else None, dir=d)
    write_units(df, paths["units"])
    write_hits(h, motifs, paths["hits"])
    motifs.to_csv(paths["families"], sep="\t", index=False)
    if write_pol:
        write_polarized(df, paths["polarized"])
    return paths, df, motifs, h


def swap_bp_tp(units_path, out_path):
    """Relabel BP <-> TP at every peak unit (bp_mask -> bp_mask ^ (2^k - 1)); hits are unchanged."""
    u = pd.read_csv(units_path, sep="\t", dtype=str, keep_default_na=False)
    peak = u["source"] == "peak"
    k = u.loc[peak, "k"].astype(int)
    u.loc[peak, "bp_mask"] = (u.loc[peak, "bp_mask"].astype(int) ^ (2 ** k - 1)).astype(str)
    u.to_csv(out_path, sep="\t", index=False)
