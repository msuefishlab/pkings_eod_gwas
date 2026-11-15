#!/usr/bin/env python3
"""tasks_ext.py -- task lists for the Stage 2b rounds, mirroring the drift test's run_fits.sh (fit / refine / nested):
  round1        4 seeds (1-4) x 5 random starts per pair x model, maxtime 2400 s
  refine R      R = 2 (fold 0.5, seeds 201-204) or 3 (fold 0.25, seeds 301-304), from the best start so far
  nested P      nested-start guard: P = 1 (seeds 351-354, after refine 2) or 2 (361-364, after refine 3); a model whose
                best LL lies below the best of its nested sub-models (ext_models.SUBMODELS, same folding) gets 4 x 5
                starts from that sub-model's best fit embedded in it, fold 0.25
Every line also carries --bounds and --theta_ref for fit_ext.py:
  round1: theta of the pair's best M0 (unfolded) / M0f (folded) start; refine: theta of the start refined;
  nested: theta of the embedded sub-model fit (the embedded spectrum is the sub-model's, so its theta is exact).
--refdir (wide-bound run): directory of the primary-bound fits; the same model at primary bounds is then an extra
nested sub-model (primary box c wide box), embedded unchanged.
Prints task lines to --out and a record of the starts to --record.
"""
import argparse
import csv
import glob
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ext_models  # noqa: E402

FIXED = {"pair", "model", "target", "seed", "start", "LL", "theta", "seconds"}


def rows_of(fitsdir, pair, model):
    return [r for f in glob.glob(os.path.join(fitsdir, f"{pair}.{model}.data.seed*.tsv"))
            for r in csv.DictReader(open(f), delimiter="\t")]


def best(fitsdir, pair, model):
    rs = rows_of(fitsdir, pair, model)
    return max(rs, key=lambda r: float(r["LL"])) if rs else None


def params_of(r):
    return {k: float(v) for k, v in r.items() if k not in FIXED}


def embed(sub, p, sup, theta, L, mode):
    """Parameters of model `sup` that reproduce sub-model `sub` at its fit p (dict); base names without 'f'."""
    lo, hi = ext_models.tc_box(theta, L, mode)
    s, S = sub.rstrip("f"), sup.rstrip("f")
    if s in ("M0", "M1", "M2"):
        T = p["T"]
        T_c = max(lo, min(math.sqrt(lo * hi), T / 2))
        q = {"nu_BP": p["nu_BP"], "nu_TP": p["nu_TP"], "T_pre": max(T - T_c, 1e-5), "T_c": T_c, "r_BP": 1.0}
        if s == "M0":
            mig = (0.0, 0.0, 0.5)
        elif s == "M1":
            mig = (p["m"], p["m"], 0.0)
        else:
            mig = (p["m_BP_from_TP"], p["m_TP_from_BP"], 0.0)
    else:
        q = {k: p[k] for k in ext_models.PNAMES[s]}
        mig = (q.pop("m_BP_from_TP", 0.0), q.pop("m_TP_from_BP", 0.0), q.pop("f_mig", 0.5))
    q.setdefault("r_TP", 1.0)
    q["m_BP_from_TP"], q["m_TP_from_BP"], q["f_mig"] = mig
    out = [q[k] for k in ext_models.PNAMES[S]]
    if "p_misid" in p:
        out.append(p["p_misid"])
    return out


def fmt(v):
    return f"{v:.8g}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["round1", "refine", "nested"])
    ap.add_argument("--n", type=int, default=0, help="refine round (2|3) or nested pass (1|2)")
    ap.add_argument("--fits", required=True)
    ap.add_argument("--sfs", required=True, help="sfs dir (for L_mask)")
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--models", required=True)
    ap.add_argument("--bounds", default="primary")
    ap.add_argument("--refdir", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--record", required=True)
    a = ap.parse_args()
    import pickle
    lines, rec = [], []
    for pair in a.pairs.split(","):
        pk = pickle.load(open(os.path.join(a.sfs, f"{pair}.pkl"), "rb"))
        for model in a.models.split(","):
            folded = model.endswith("f")
            f = "f" if folded else ""
            L = float(pk["L_mask_alt"] if folded else pk["L_mask"])
            if a.mode == "round1":
                b0 = best(a.fits, pair, "M0" + f)
                if b0 is None:
                    sys.exit(f"no M0{f} fit for {pair}: theta_ref undefined")
                th = float(b0["theta"])
                for s in (1, 2, 3, 4):
                    lines.append(f"--pair {pair} --model {model} --starts 5 --seed {s} --maxtime 2400 "
                                 f"--bounds {a.bounds} --theta_ref {th:.8g}")
                rec.append(f"{pair}\t{model}\ttheta_ref=M0{f}\t{th:.8g}")
            elif a.mode == "refine":
                fold = {2: 0.5, 3: 0.25}[a.n]
                b = best(a.fits, pair, model)
                if b is None:
                    continue
                p0 = ",".join(b[k] for k in b if k not in FIXED)
                for s in (1, 2, 3, 4):
                    lines.append(f"--pair {pair} --model {model} --starts 5 --seed {a.n * 100 + s} --p0 {p0} "
                                 f"--fold {fold} --maxtime 2400 --bounds {a.bounds} --theta_ref {b['theta']}")
                rec.append(f"{pair}\t{model}\t{p0}\tLL={b['LL']}\ttheta={b['theta']}")
            else:
                base = {1: 350, 2: 360}[a.n]
                B = best(a.fits, pair, model)
                cands = []
                for sub in ext_models.SUBMODELS[model.rstrip("f")]:
                    bs = best(a.fits, pair, sub + f)
                    if bs is not None:
                        cands.append((sub + f, bs, "same-run"))
                if a.refdir:
                    bs = best(a.refdir, pair, model)
                    if bs is not None:
                        cands.append((model, bs, "primary-bounds"))
                if not cands or B is None:
                    continue
                sub, bs, src = max(cands, key=lambda c: float(c[1]["LL"]))
                if float(B["LL"]) >= float(bs["LL"]):
                    continue
                th = float(bs["theta"])
                if src == "primary-bounds":
                    p0 = [float(bs[k]) for k in bs if k not in FIXED]
                else:
                    p0 = embed(sub, params_of(bs), model, th, L, a.bounds)
                p0s = ",".join(fmt(v) for v in p0)
                for s in (1, 2, 3, 4):
                    lines.append(f"--pair {pair} --model {model} --starts 5 --seed {base + s} --p0 {p0s} "
                                 f"--fold 0.25 --maxtime 2400 --bounds {a.bounds} --theta_ref {th:.8g}")
                rec.append(f"{pair}\t{model}\tbelow {sub} ({src}) by {float(bs['LL']) - float(B['LL']):.1f}\t{p0s}")
    with open(a.out, "w") as o:
        o.write("".join(l + "\n" for l in lines))
    with open(a.record, "w") as o:
        o.write("".join(l + "\n" for l in rec))
    print(f"{a.mode} {a.n}: {len(lines)} tasks")


if __name__ == "__main__":
    main()
