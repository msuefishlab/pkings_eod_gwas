#!/usr/bin/env python3
"""gen_tasks.py -- Stage 3 task lists (lines '<out_check> <tab> <kind> <tab> <args>' for task.sh).

  families <models,comma>          step 3: per model and simulation/23e grid value (i_ratio 1..12):
        cal  100 chunks = 100 calibration blocks drawn once, default_rng((20261109, 1)).choice(468, 100, replace=False)
             (g = 1..100, chunk c on the c-th drawn block)  -- reference curves for the calibration fit (criterion P1)
        ins  100 chunks = 100 M-sum blocks drawn once, default_rng((20261109, 2)).choice(860, 100, replace=False)
             (g = 101..200)  -- reference curves for the in-sample continuity form (simulation/23e)
        ho   held-out scoring family: 5 replicate sets (g = 1001..1005) x all 392 held-out blocks (split in 4 parts)
  genomes <model> <i_ratio> <ratio> <g_from> <g_to>   whole M-sum genomes (25 chromosome tasks each) + post lines
Usage: gen_tasks.py <mode> ... <out_tasks.txt> [<out_post.txt>]
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim_lib as S  # noqa: E402

O = os.environ["O"]
BL = S.load_blocks(os.environ["S3_MASK"], os.environ["S3_SPLIT"])
CAL = [b["j"] for b in BL if b["set"] == "calibration"]
HO = [b["j"] for b in BL if b["set"] == "heldout"]
ALL = [b["j"] for b in BL]
assert len(CAL) == 468 and len(HO) == 392 and len(ALL) == 860
mode = sys.argv[1]
lines, post = [], []


def draws():
    cal = [CAL[i] for i in np.random.default_rng((20261109, 1)).choice(len(CAL), 100, replace=False)]
    ins = [ALL[i] for i in np.random.default_rng((20261109, 2)).choice(len(ALL), 100, replace=False)]
    return cal, ins


if mode == "families":
    models = sys.argv[2].split(","); out = sys.argv[3]
    SPLIT_CAL = os.environ.get("SPLIT_CAL", "0") == "1"   # load balancing only (alternatives run); chunks unchanged
    cal, ins = draws()
    for k in range(12, 0, -1):                         # high r/mu first (longest tasks), across all models
        for m in models:
            d = f"{O}/fam/{m}"; os.makedirs(f"{d}/lists", exist_ok=True)
            r = S.GRID[k - 1]
            npart = 4 if SPLIT_CAL else 1
            units = {}
            for p, part in enumerate(np.array_split(np.arange(100), npart)):
                sfx = f"_p{p}" if npart > 1 else ""
                units[f"cal{sfx}"] = [(c + 1, cal[c]) for c in part]
                units[f"ins{sfx}"] = [(101 + c, ins[c]) for c in part]
            for s_ in range(1, 6):
                for p, part in enumerate(np.array_split(np.array(HO), 4)):
                    units[f"ho_s{s_}_p{p}"] = [(1000 + s_, int(j)) for j in part]
            for fam, us in units.items():
                lst = f"{d}/lists/{fam}_k{k}.tsv"
                with open(lst, "w") as fo:
                    fo.write("".join(f"{g}\t{j}\n" for g, j in us))
                o = f"{d}/{fam}_k{k}.npz"
                lines.append(f"{o}\tsim\tchunks {m} {k} {r:g} {lst} {o}")
    with open(f"{O}/fam/draws.tsv", "w") as fo:
        fo.write("family\tchunk\tg\tj\tchrom\tblock\tset\n")
        for c, j in enumerate(cal):
            fo.write(f"cal\t{c + 1}\t{c + 1}\t{j}\t{BL[j]['chrom']}\t{BL[j]['block']}\t{BL[j]['set']}\n")
        for c, j in enumerate(ins):
            fo.write(f"ins\t{c + 1}\t{101 + c}\t{j}\t{BL[j]['chrom']}\t{BL[j]['block']}\t{BL[j]['set']}\n")
elif mode == "genomes":
    m, ir, ratio, g0, g1, out, outp = sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), \
        int(sys.argv[6]), sys.argv[7], sys.argv[8]
    # chromosomes in decreasing block count (longest first)
    nb = {}
    for b in BL:
        nb[b["chrom"]] = nb.get(b["chrom"], 0) + 1
    chs = sorted(nb, key=lambda c: -nb[c])
    for g in range(g0, g1 + 1):
        gd = f"{O}/genomes/{m}/r{ir}/g{g}"
        for ch in chs:
            lines.append(f"{gd}/seg/inter.{ch}.npz\tsim\tgenome {m} {ir} {ratio!r} {g} {ch} {gd}/seg")
        post.append(f"{gd}/chromld.npz\tpost\t{gd} {os.environ['S3_SPLIT']} {O}/samples_def_BAM_TP1.tsv "
                    f"{O}/pairs_def_BAM_TP1.tsv")
    with open(outp, "a") as fo:
        fo.write("".join(x + "\n" for x in post))
    with open(out, "a") as fo:
        fo.write("".join(x + "\n" for x in lines))
    print(len(lines), "genome tasks,", len(post), "post tasks")
    sys.exit(0)
else:
    sys.exit("mode?")
with open(out, "a") as fo:
    fo.write("".join(x + "\n" for x in lines))
print(len(lines), "tasks")
