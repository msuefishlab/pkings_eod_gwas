#!/usr/bin/env python3
"""family_summary.py -- pool the step-3 chunk files into compact sums.
For model m, family fam (cal | ins | ho_s1..ho_s5), grid index k, treatment t, population pop:
  <m>|<fam>|<k>|<t>|<pop>|cnt, ...|tot   summed n_pairs and sum r^2 per bin over chunks (the 23e pooling)
  <m>|<fam>|<k>|<t>|<pop>|se             MC SE of the pooled curve, bootstrap over chunks (200 resamples, rng 20261113)
  <m>|<fam>|<k>|time                     [cpu s simulation, cpu s scoring, n chunks]
Usage: family_summary.py <fam_dir> <models,comma> <out.npz>
"""
import glob, os, re, sys
import numpy as np

fd, models, out = sys.argv[1], sys.argv[2].split(","), sys.argv[3]
acc, raw, tim = {}, {}, {}
for m in models:
    for fn in sorted(glob.glob(f"{fd}/{m}/*_k*.npz")):
        if ".tmp." in fn:
            continue
        b = os.path.basename(fn)[:-4]
        mo = re.match(r"(cal|ins|ho_s(\d))(?:_p\d)?_k(\d+)$", b)
        fam, k = (mo.group(1), int(mo.group(3)))
        z = np.load(fn)
        t = tim.setdefault(f"{m}|{fam}|{k}|time", np.zeros(3))
        t += [z["t_sim"].sum(), z["t_score"].sum(), len(z["g"])]
        for key in z.files:
            if key.endswith("|cnt") or key.endswith("|tot"):
                base = key.rsplit("|", 1)[0]
                raw.setdefault((m, fam, k, base, key[-3:]), []).append(z[key])
rng = np.random.default_rng(20261113)
for (m, fam, k, base, kind), L in raw.items():
    A = np.vstack(L)
    acc[f"{m}|{fam}|{k}|{base}|{kind}"] = A.sum(0)
    if kind == "cnt":
        T = np.vstack(raw[(m, fam, k, base, "tot")])
        bs = []
        for _ in range(200):
            i = rng.integers(0, A.shape[0], A.shape[0])
            with np.errstate(invalid="ignore", divide="ignore"):
                bs.append(T[i].sum(0) / A[i].sum(0))
        acc[f"{m}|{fam}|{k}|{base}|se"] = np.nanstd(bs, axis=0)
acc.update(tim)
np.savez_compressed(out + ".tmp.npz", **acc)
os.replace(out + ".tmp.npz", out)
print("keys", len(acc), "chunk-files", sum(int(v[2] > 0) for v in tim.values()))
