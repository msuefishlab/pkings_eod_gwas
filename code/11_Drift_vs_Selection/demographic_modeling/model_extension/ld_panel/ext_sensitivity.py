#!/usr/bin/env python3
"""ext_sensitivity.py -- shift of every panel statistic when the admixed 2019 Apassa wobbles or
close relatives are removed.

For each contrast (variant vs reference, same statistic and chromosome set):
  shift = variant - reference; shift_halfbands = shift / ((hi_ref - lo_ref) / 2) (the reference's chromosome-bootstrap
  half-width); outside_ref_band = variant value outside [lo_ref, hi_ref];
  control_shift_min/max/mean = the same shift for the size-matched random-removal controls (same number of fish removed,
  drawn with a fixed seed from unflagged fish of the same locality); net_shift = shift - mean control shift (the
  effect of removing THESE fish rather than any fish of that locality, i.e. net of sample size and of locality mix);
  net_shift_halfbands = net_shift / reference half-width; outside_control_range = shift outside [min, max] of the
  control shifts (5 controls: the variant is more extreme than every control).
Usage: ext_sensitivity.py <panel.tsv> [<panel.tsv> ...] -o <out.tsv>
  (first file wins for duplicated (statistic, sample, set) rows; pass the corrected admix2 panel first)
"""
import sys
import numpy as np
import pandas as pd

args = sys.argv[1:]
out = args[args.index("-o") + 1]
fns = [a for a in args[:args.index("-o")]]
P = pd.concat([pd.read_csv(f, sep="\t") for f in fns]).drop_duplicates(["statistic", "sample", "set"])
C = [  # contrast, reference unit, variant unit, control units
    # (a) admixture, PRIMARY: TP1 minus all ten 2019 Apassa wobbles (metadata POP=APA, Date 2019)
    ("a_admixture_TP1_noAdmix", "TP1", "TP1_noAdmix", [f"TP1_minus10A09_r{k}" for k in range(1, 6)]),
    ("a_admixture_TP1_noAdmix", "BAM|TP1", "BAM|TP1_noAdmix", [f"BAM|TP1_minus10A09_r{k}" for k in range(1, 6)]),
    # (a2) sensitivity: TP1 minus only the seven 2019-017 fish
    ("a2_admixture_TP1_no017", "TP1", "TP1_no017", [f"TP1_minus7A09_r{k}" for k in range(1, 6)]),
    ("a2_admixture_TP1_no017", "BAM|TP1", "BAM|TP1_no017", [f"BAM|TP1_minus7A09_r{k}" for k in range(1, 6)]),
    # (b) kinship: drop BAM_175 (second-degree pair BAM_167/BAM_175); third-degree pairs kept
    ("b_kinship2_BAM_noKin", "BAM", "BAM_noKin", [f"BAM_minus1_r{k}" for k in range(1, 6)]),
    ("b_kinship2_BAM_noKin", "BAM|TP1", "BAM_noKin|TP1", [f"BAM_minus1_r{k}|TP1" for k in range(1, 6)]),
    # extra, beyond the admixture and kinship rules above: third-degree removal (BEN_726, Bengue)
    ("extra_kinship3_TP1_noKin3", "TP1", "TP1_noKin3", [f"TP1_minus1BEN_r{k}" for k in range(1, 6)]),
    ("extra_kinship3_TP1_noKin3", "BAM|TP1", "BAM|TP1_noKin3", [f"BAM|TP1_minus1BEN_r{k}" for k in range(1, 6)]),
]
ix = P.set_index(["statistic", "sample", "set"])
rows = []
for name, ref, var, ctrls in C:
    for (stat, s, st), r in ix.loc[(slice(None), ref, slice(None)), :].iterrows():
        if (stat, var, st) not in ix.index:
            continue
        v = ix.loc[(stat, var, st)]
        hw = (r["hi"] - r["lo"]) / 2
        sh = v["value"] - r["value"]
        cs = np.array([ix.loc[(stat, c, st)]["value"] - r["value"] for c in ctrls if (stat, c, st) in ix.index])
        rows.append(dict(contrast=name, statistic=stat, set=st, reference=ref, variant=var,
                         ref_value=r["value"], ref_lo=r["lo"], ref_hi=r["hi"], variant_value=v["value"],
                         variant_lo=v["lo"], variant_hi=v["hi"], shift=sh,
                         shift_halfbands=sh / hw if hw > 0 else np.nan,
                         outside_ref_band=int(not (r["lo"] <= v["value"] <= r["hi"])),
                         n_controls=len(cs),
                         control_shift_min=cs.min() if len(cs) else np.nan,
                         control_shift_max=cs.max() if len(cs) else np.nan,
                         control_shift_mean=cs.mean() if len(cs) else np.nan,
                         net_shift=sh - cs.mean() if len(cs) else np.nan,
                         net_shift_halfbands=(sh - cs.mean()) / hw if len(cs) and hw > 0 else np.nan,
                         outside_control_range=int(sh < cs.min() or sh > cs.max()) if len(cs) else -1))
pd.DataFrame(rows).to_csv(out, sep="\t", index=False, float_format="%.6g")
print(len(rows), "rows")
