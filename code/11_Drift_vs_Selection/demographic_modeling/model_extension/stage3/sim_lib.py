"""sim_lib.py -- model extension, Stage 3: physical-coordinate simulation of one M-sum block and its scoring.

Simulation (as coord_pilot arm b''): the 1-Mb block plus 100-kb flanks on each side (1.2 Mb), uniform
recombination r = (physical r/mu) x mu per bp, mu = 1e-8, msprime sim_ancestry (67 BP + 47 TP diploids, ploidy 2) and
sim_mutations with the BinaryMutationModel on the discrete genome; physical position = block start - 100 kb + site + 1.
Kept: mutations at every position of the block that has a record in the all-sites VCF (primary LD site set); each kept
site carries that position's observed per-individual codes (records_full.py): missingness, DP bin, strict-genotype flag; callable
flag of the position. Flanks and positions without a record are discarded.

Demography: baseline models -- the 23b/04c function, verbatim (params.as_run.agree.json, mu = 1e-8 arm);
alternatives -- piecewise-constant epochs from the 01_fitted_history.py conversion (build_alt).

Treatments (all on the SAME simulated genotypes):
  primary   standard missingness; LD on all kept sites (primary LD site set); pi / dxy / Fst / Delta on callable kept sites
  callable  standard missingness; LD on callable kept sites only (sensitivity)
  gqdrop    each called heterozygote of individual i becomes REF or ALT homozygote (equal odds) with prob d(i, b),
            b = observed DP bin of i at that position (uncalibrated sensitivity); standard missingness
  d019      strict-genotype missingness (DP < 8 or GQ/RGQ < 20 or missing)
  eps0.005, eps0.01  symmetric genotype error: each genotype with prob eps replaced by one of the other two genotypes
            (equal odds) (no empirical basis); standard missingness
Polarization: unfolded models flip each kept site's derived/ancestral label with p_misid (pol = -1);
folded models pol = +1 everywhere and polarized statistics are not used.
Seeds: one genome / chunk seed s (sim_index.tsv rule); per block j (global mask index) SeedSequence([s, j]) gives the
ancestry and mutation seeds; treatment rngs default_rng([s, j, k]) with k = 1 misid, 2 gqdrop, 3 eps0.005, 4 eps0.01.
"""
import csv, json, os
import numpy as np
import msprime

MU = 1e-8
FLANK = 100_000
N_BP, N_TP = 67, 47
NI = N_BP + N_TP
SAMPLES = {"BAM": np.arange(N_BP), "TP1": np.arange(N_BP, NI)}
PAIRS = [("BAM", "TP1")]
GRID = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0]
TREATMENTS = ["primary", "callable", "gqdrop", "d019", "eps0.005", "eps0.01"]
SEED0 = 20261109


def genome_seed(i_model, i_ratio, g):
    return SEED0 + 1_000_000 * int(i_model) + 10_000 * int(i_ratio) + int(g)


def block_seeds(s, j):
    a, m = np.random.SeedSequence([int(s), int(j)]).generate_state(2)
    return int(a) % (2 ** 31 - 2) + 1, int(m) % (2 ** 31 - 2) + 1


def load_blocks(mask_fn, split_fn):
    split = {r["chrom"]: r["set"] for r in csv.DictReader(open(split_fn), delimiter="\t")}
    out = []
    for j, r in enumerate(csv.DictReader(open(mask_fn), delimiter="\t")):
        out.append(dict(j=j, chrom=r["chrom"], block=int(r["block"]), f=float(r["f"]), L=int(r["L_callable"]),
                        set=split[r["chrom"]]))
    return out


def load_models(index_fn):
    """model rows of sim_index.tsv; alternatives get p_misid / phys_json from alt_params.tsv (env S3_ALT) when present"""
    M = {r["model"]: r for r in csv.DictReader(open(index_fn), delimiter="\t") if r.get("kind") == "model"}
    alt = os.environ.get("S3_ALT", "")
    if alt and os.path.exists(alt):
        for r in csv.DictReader(open(alt), delimiter="\t"):
            if r["model"] in M:
                M[r["model"]]["p_misid"] = r["p_misid"]
                M[r["model"]]["phys_json"] = r["phys_json"]
    return M


# ---------------------------------------------------------------- demography
def demography_baseline(P, model):
    # verbatim from 23b_sim_ld_decay.py (as 04c_worker_msprime_agree.py)
    d = msprime.Demography()
    d.add_population(name="BP", initial_size=P["N_BP_end"])
    d.add_population(name="TP", initial_size=P["N_TP_end"])
    d.add_population(name="ANC", initial_size=P["N_anc"])
    if model.rstrip("f") == "M3":
        T = P["T_gen"]
        d.populations[0].growth_rate = np.log(P["N_BP_end"] / P["N_BP_start"]) / T
        d.populations[1].growth_rate = np.log(P["N_TP_end"] / P["N_TP_start"]) / T
    if P["m_BP_from_TP_pergen"] > 0:
        d.set_migration_rate(source="BP", dest="TP", rate=P["m_BP_from_TP_pergen"])
    if P["m_TP_from_BP_pergen"] > 0:
        d.set_migration_rate(source="TP", dest="BP", rate=P["m_TP_from_BP_pergen"])
    d.add_population_split(time=P["T_gen"], derived=["BP", "TP"], ancestral="ANC")
    return d


def demography_alt(P):
    """Alternatives CB / CBT / CBTM (+f), from physical parameters (mu = 1e-8 arm,
    01_fitted_history.py conversion as in stage2b ext_models.physical): present sizes N_BP_2 / N_TP_2; at T_c_gen before
    present the sizes become N_BP_1 / N_TP_1; migration (dadi m12 = into BP from TP, forwards: msprime source BP ->
    dest TP backwards, as 23b) acts from the split until T_m_gen before present, i.e. backwards on [T_m_gen, T_s_gen];
    split at T_s_gen into ANC (N_anc)."""
    d = msprime.Demography()
    d.add_population(name="BP", initial_size=P["N_BP_2"])
    d.add_population(name="TP", initial_size=P["N_TP_2"])
    d.add_population(name="ANC", initial_size=P["N_anc"])
    mb, mt = float(P.get("m_BP_from_TP_pergen", 0) or 0), float(P.get("m_TP_from_BP_pergen", 0) or 0)
    tm_raw = P.get("T_m_gen", None)
    Tm = float(tm_raw) if (tm_raw not in (None, "") and (mb > 0 or mt > 0)) else float(P["T_s_gen"])
    if (mb > 0 or mt > 0) and Tm <= 0:
        if mb > 0:
            d.set_migration_rate(source="BP", dest="TP", rate=mb)
        if mt > 0:
            d.set_migration_rate(source="TP", dest="BP", rate=mt)
    elif (mb > 0 or mt > 0) and Tm < P["T_s_gen"]:
        if mb > 0:
            d.add_migration_rate_change(time=Tm, rate=mb, source="BP", dest="TP")
        if mt > 0:
            d.add_migration_rate_change(time=Tm, rate=mt, source="TP", dest="BP")
    if P["T_c_gen"] < P["T_s_gen"]:
        d.add_population_parameters_change(time=P["T_c_gen"], initial_size=P["N_BP_1"], population="BP")
        d.add_population_parameters_change(time=P["T_c_gen"], initial_size=P["N_TP_1"], population="TP")
    d.add_population_split(time=P["T_s_gen"], derived=["BP", "TP"], ancestral="ANC")
    d.sort_events()
    return d


def model_demography(spec, params_json):
    if spec["source"] == "baseline":
        P = json.load(open(params_json))[f"{spec['model']}|1e-08"]
        return demography_baseline(P, spec["model"])
    P = json.loads(spec["phys_json"])
    return demography_alt(P)


# ---------------------------------------------------------------- records cache
class Records:
    """per-block record arrays from a cache dir (rpos.npy, code.npy, cal.npy per block; mmap)."""

    def __init__(self, cache):
        self.cache = cache

    def get(self, ch, b):
        p = os.path.join(self.cache, f"{ch}_{b}")
        return (np.load(p + ".rpos.npy"), np.load(p + ".code.npy", mmap_mode="r"), np.load(p + ".cal.npy"))


def load_dropout(fn):
    D = np.zeros((NI, 5))
    for r in csv.DictReader(open(fn), delimiter="\t"):
        D[int(r["i"]), int(r["bin"])] = float(r["d"])
    return D


# ---------------------------------------------------------------- simulate
def simulate_block(dem, ratio, blk, rec, sa, sm):
    """returns X (n_kept, 114) int8 true dosages, pos (n_kept,) physical, code rows (n_kept, 114) uint8, cal (n_kept,)"""
    S0 = blk["block"] * 1_000_000
    ts = msprime.sim_ancestry(samples={"BP": N_BP, "TP": N_TP}, demography=dem, sequence_length=1_000_000 + 2 * FLANK,
                              recombination_rate=ratio * MU, ploidy=2, random_seed=sa)
    ts = msprime.sim_mutations(ts, rate=MU, model=msprime.BinaryMutationModel(), discrete_genome=True,
                               random_seed=sm)
    H = ts.genotype_matrix()
    phys = S0 - FLANK + ts.tables.sites.position.astype(np.int64) + 1
    rpos, code, cal = rec
    j = np.clip(np.searchsorted(rpos, phys), 0, max(len(rpos) - 1, 0))
    keep = (rpos[j] == phys) if len(rpos) else np.zeros(len(phys), bool)
    H = H[keep]
    X = (H[:, 0::2] + H[:, 1::2]).astype(np.int8)
    jj = j[keep]
    return X, phys[keep], np.asarray(code[jj]), cal[jj], ts.num_trees


def apply_treatment(t, X, code, rng=None, D=None):
    miss = ((code >> 3) & 1).astype(bool)
    if t in ("primary", "callable"):
        return np.where(miss, -1, X).astype(np.int8)
    if t == "d019":
        return np.where(((code >> 4) & 1).astype(bool), -1, X).astype(np.int8)
    if t == "gqdrop":
        p = D[np.arange(NI)[None, :], (code & 7).astype(int)]
        u = rng.random(X.shape); v = rng.random(X.shape)
        Y = np.where((X == 1) & (u < p), np.where(v < 0.5, 0, 2), X)
        return np.where(miss, -1, Y).astype(np.int8)
    if t.startswith("eps"):
        eps = float(t[3:])
        u = rng.random(X.shape); v = rng.random(X.shape)
        Y = np.where(u < eps, (X + 1 + (v < 0.5)) % 3, X)
        return np.where(miss, -1, Y).astype(np.int8)
    raise ValueError(t)


TK = {"gqdrop": 2, "eps0.005": 3, "eps0.01": 4}


def ld_curve(L, P, G, pos, f, idx):
    keep, Xs, Ms = P.ld_keep(G[:, idx])
    if keep.sum() >= 2:
        cnt, tot = L.binned_r2(pos[keep], Xs[keep], Ms[keep], f)
    else:
        cnt, tot = np.zeros(L.NBIN, np.int64), np.zeros(L.NBIN)
    return np.asarray(cnt, float), np.asarray(tot, float), int(keep.sum())
