#!/usr/bin/env python3
"""sim_task.py -- one Stage 3 simulation task (see sim_lib.py for the simulation and the treatments).

Modes
  genome <model> <i_ratio> <ratio> <g> <chrom> <out_dir>
      one chromosome of whole M-sum genome replicate g: every M-sum block of <chrom> in mask order, primary treatment,
      panel_lib.segment_stats (pi, dxy, Fst, Delta ladder, LD BAM/TP1) per block -> <out_dir>/seg.<chrom>.tsv and the
      interchromosomal draw (panel_lib.draw_interchrom, Stage 1 seed rule) -> <out_dir>/inter.<chrom>.npz, in exactly the
      format of ld_panel/ext_panel_chrom.py, so ld_panel/ext_aggregate.py runs on the genome directory unchanged.
  chunks <model> <i_ratio> <ratio> <chunk_list.tsv> <out.npz>
      reference / scoring-family chunks: lines "g <tab> j"; one block-replicate per line, LD only:
      BAM under all six treatments, TP1 under primary.
  pilot <i_block_pilot> <i_ratio_pilot> <pilot_blocks.tsv> <p_misid> <out.npz>
      reproduction of coord_pilot arm b'' (pilot_sim_a2.py bpp, 100 replicates, pilot seeds, M2) through this code:
      LD of bpp (all record sites) and b_chk (callable subset), BAM and TP1.
Env: S3_MASK, S3_SPLIT, S3_INDEX, S3_PARAMS, S3_CACHE, S3_DROPOUT, PANEL_DIR, LD23_DIR
"""
import csv, json, os, sys, time
import numpy as np

sys.path.insert(0, os.environ["PANEL_DIR"])
import panel_lib as P   # noqa: E402  (unchanged; imports ld23_lib unchanged)
import ld23_lib as L    # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim_lib as S     # noqa: E402

mode = sys.argv[1]
BLOCKS = S.load_blocks(os.environ["S3_MASK"], os.environ["S3_SPLIT"])
REC = S.Records(os.environ["S3_CACHE"])


def write_seg(fo, ch, b, st):
    for (kind, name), v in st.items():
        if np.ndim(v):
            for i, x in enumerate(v):
                fo.write(f"{ch}\t{b}\t{kind}\t{name}\t{i}\t{x:.10g}\n")
        else:
            fo.write(f"{ch}\t{b}\t{kind}\t{name}\t\t{v:.10g}\n")


if mode == "genome":
    model, i_ratio, ratio, g, chrom, out = sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), \
        sys.argv[6], sys.argv[7]
    spec = S.load_models(os.environ["S3_INDEX"])[model]
    dem = S.model_demography(spec, os.environ["S3_PARAMS"])
    pmis = float(spec["p_misid"])
    s = S.genome_seed(spec["i_model"], i_ratio, g)
    os.makedirs(out, exist_ok=True)
    rows = {n: ([], []) for n in S.SAMPLES}
    tmp = os.path.join(out, f"seg.{chrom}.tsv.tmp")
    t0 = time.process_time(); nb = 0
    with open(tmp, "w") as fo:
        for blk in (x for x in BLOCKS if x["chrom"] == chrom):
            sa, sm = S.block_seeds(s, blk["j"])
            X, pos, code, cal, _ = S.simulate_block(dem, ratio, blk, REC.get(chrom, blk["block"]), sa, sm)
            r1 = np.random.default_rng([s, blk["j"], 1])
            pol = np.where(r1.random(len(pos)) < pmis, -1, 1).astype(np.int8)
            G = S.apply_treatment("primary", X, code)
            st = P.segment_stats(G, pos, blk["f"], blk["L"], cal, pol, S.SAMPLES, S.PAIRS)
            write_seg(fo, chrom, blk["block"], st)
            for n, ix in S.SAMPLES.items():
                keep, Xk, Mk = P.ld_keep(G[:, ix])
                rows[n][0].append(Xk[keep]); rows[n][1].append(Mk[keep])
            nb += 1
    os.replace(tmp, os.path.join(out, f"seg.{chrom}.tsv"))
    c = int(chrom[3:])
    arr = {}
    for n, (xs, ms) in rows.items():
        Xc = np.concatenate(xs) if xs else np.zeros((0, len(S.SAMPLES[n])), np.int8)
        Mc = np.concatenate(ms) if ms else np.zeros((0, len(S.SAMPLES[n])), np.int8)
        d = P.draw_interchrom({c: (Xc, Mc)})[c]
        arr[f"{n}__X"], arr[f"{n}__M"], arr[f"{n}__Xp"], arr[f"{n}__Mp"] = (a.astype(np.int8) for a in d)
        arr[f"{n}__m"] = np.array(Xc.shape[0])
    np.savez_compressed(os.path.join(out, f"inter.{chrom}.tmp.npz"), **arr)
    os.replace(os.path.join(out, f"inter.{chrom}.tmp.npz"), os.path.join(out, f"inter.{chrom}.npz"))
    print(f"genome {model} r{i_ratio} ratio={ratio:g} g={g} {chrom} blocks={nb} cpu={time.process_time() - t0:.1f}s",
          flush=True)

elif mode == "chunks":
    model, i_ratio, ratio, lst, out = sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), sys.argv[5], sys.argv[6]
    spec = S.load_models(os.environ["S3_INDEX"])[model]
    dem = S.model_demography(spec, os.environ["S3_PARAMS"])
    D = S.load_dropout(os.environ["S3_DROPOUT"])
    units = [tuple(int(x) for x in l.split()) for l in open(lst) if l.strip()]
    KEYS = [(t, "BAM") for t in S.TREATMENTS] + [("primary", "TP1")]
    acc = {k: ([], [], []) for k in KEYS}
    tsim, tsc, t0 = [], [], time.process_time()
    for g, j in units:
        blk = BLOCKS[j]
        s = S.genome_seed(spec["i_model"], i_ratio, g)
        sa, sm = S.block_seeds(s, j)
        ta = time.process_time()
        X, pos, code, cal, _ = S.simulate_block(dem, ratio, blk, REC.get(blk["chrom"], blk["block"]), sa, sm)
        tb = time.process_time()
        for t in S.TREATMENTS:
            rng = np.random.default_rng([s, j, S.TK[t]]) if t in S.TK else None
            G = S.apply_treatment(t, X, code, rng, D)
            Gt, pt = (G[cal], pos[cal]) if t == "callable" else (G, pos)
            pops = ["BAM", "TP1"] if t == "primary" else ["BAM"]
            for pop in pops:
                cnt, tot, nk = S.ld_curve(L, P, Gt, pt, blk["f"], S.SAMPLES[pop])
                a = acc[(t, pop)]; a[0].append(cnt); a[1].append(tot); a[2].append(nk)
        tsim.append(tb - ta); tsc.append(time.process_time() - tb)
    arrs = {"g": np.array([u[0] for u in units]), "j": np.array([u[1] for u in units]),
            "t_sim": np.array(tsim), "t_score": np.array(tsc)}
    for (t, pop), (c_, s_, n_) in acc.items():
        arrs[f"{t}|{pop}|cnt"] = np.vstack(c_) if c_ else np.zeros((0, L.NBIN))
        arrs[f"{t}|{pop}|tot"] = np.vstack(s_) if s_ else np.zeros((0, L.NBIN))
        arrs[f"{t}|{pop}|nsnp"] = np.array(n_)
    meta = dict(model=model, i_ratio=i_ratio, ratio=ratio, n=len(units))
    tmp = out[:-4] + ".tmp.npz"
    np.savez_compressed(tmp, meta=json.dumps(meta), **arrs)
    os.replace(tmp, out)
    print(f"chunks {model} r{i_ratio} ratio={ratio:g} n={len(units)} cpu={time.process_time() - t0:.1f}s "
          f"sim={np.sum(tsim):.1f}s score={np.sum(tsc):.1f}s", flush=True)

elif mode == "pilot":
    ib, ir, pb, pmis, out = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], float(sys.argv[5]), sys.argv[6]
    spec = S.load_models(os.environ["S3_INDEX"])["M2"]
    dem = S.model_demography(spec, os.environ["S3_PARAMS"])
    pr = [r for r in csv.DictReader(open(pb), delimiter="\t")][ib]
    blk = next(x for x in BLOCKS if x["chrom"] == pr["chrom"] and x["block"] == int(pr["block"]))
    assert abs(blk["f"] - float(pr["f"])) < 1e-9 and blk["L"] == int(pr["L_callable"])
    ratio = [0.25, 0.5, 1.0, 1.4, 2.5][ir]
    rec = REC.get(blk["chrom"], blk["block"])
    acc = {}
    for rep in range(1, 101):
        seed = 20261110 + 1_000_000 * 2 + 100_000 * ir + 1_000 * ib + rep
        X, pos, code, cal, _ = S.simulate_block(dem, ratio, blk, rec, seed, seed + 50_000_000)
        G = S.apply_treatment("primary", X, code)
        for v, (Gt, pt) in (("bpp", (G, pos)), ("b_chk", (G[cal], pos[cal]))):
            for pop in ("BAM", "TP1"):
                cnt, tot, nk = S.ld_curve(L, P, Gt, pt, blk["f"], S.SAMPLES[pop])
                acc.setdefault(f"{v}|ld_cnt|{pop}", []).append(cnt)
                acc.setdefault(f"{v}|ld_tot|{pop}", []).append(tot)
    tmp = out[:-4] + ".tmp.npz"
    np.savez_compressed(tmp, **{k: np.vstack(v) for k, v in acc.items()})
    os.replace(tmp, out)
    print(f"pilot ib={ib} ir={ir} {blk['chrom']}:{blk['block']} done", flush=True)
elif mode == "sfs":
    # alternative-model unit check: joint SFS at callable positions, drift-test projection (60 x 60
    # hypergeometric, derived = simulated ALT, flipped with p_misid for unfolded models), per-position missingness.
    from math import comb
    model, i_ratio, ratio, lst, out = sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), sys.argv[5], sys.argv[6]
    spec = S.load_models(os.environ["S3_INDEX"])[model]
    dem = S.model_demography(spec, os.environ["S3_PARAMS"])
    pmis = float(spec["p_misid"] or 0)
    units = [tuple(int(x) for x in l.split()) for l in open(lst) if l.strip()]
    NPJ = 60
    fs = np.zeros((NPJ + 1, NPJ + 1)); Lsum = 0; nsites = 0; cache = {}; per_fs, per_L = [], []

    def proj(a, n):
        key = (int(a), int(n))
        if key not in cache:
            n_, a_ = int(n), int(a)          # hypergeometric projection (dadi's), exact integers
            cache[key] = np.array([comb(a_, k) * comb(n_ - a_, NPJ - k) for k in range(NPJ + 1)], float) / comb(n_, NPJ)
        return cache[key]
    for g, j in units:
        blk = BLOCKS[j]
        s = S.genome_seed(spec["i_model"], i_ratio, g)
        sa, sm = S.block_seeds(s, j)
        rec = REC.get(blk["chrom"], blk["block"])
        X, pos, code, cal, _ = S.simulate_block(dem, ratio, blk, rec, sa, sm)
        r1 = np.random.default_rng([s, j, 1])
        flip = r1.random(len(pos)) < pmis
        G = S.apply_treatment("primary", X, code)[cal]; flip = flip[cal]
        Lsum += int(rec[2].sum())
        Gb, Gt = G[:, S.SAMPLES["BAM"]], G[:, S.SAMPLES["TP1"]]
        nb_, nt_ = 2 * (Gb >= 0).sum(1), 2 * (Gt >= 0).sum(1)
        ab = np.where(Gb >= 0, Gb, 0).sum(1); at = np.where(Gt >= 0, Gt, 0).sum(1)
        ab = np.where(flip, nb_ - ab, ab); at = np.where(flip, nt_ - at, at)
        ok = (nb_ >= NPJ) & (nt_ >= NPJ) & ((ab + at) > 0)        # corners are masked at comparison
        f1 = np.zeros_like(fs)
        for x, y, u, v in zip(ab[ok], nb_[ok], at[ok], nt_[ok]):
            f1 += np.outer(proj(x, y), proj(u, v))
        fs += f1; per_fs.append(f1.astype(np.float32)); per_L.append(int(rec[2].sum()))
        nsites += int(ok.sum())
    tmp = out[:-4] + ".tmp.npz"
    np.savez_compressed(tmp, fs=fs, L=Lsum, nsites=nsites, n_units=len(units), per_fs=np.array(per_fs), per_L=np.array(per_L))
    os.replace(tmp, out)
    print(f"sfs {model} n={len(units)} L={Lsum} sites={nsites}", flush=True)
else:
    sys.exit(f"unknown mode {mode}")
