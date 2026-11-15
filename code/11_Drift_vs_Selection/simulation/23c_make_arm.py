#!/usr/bin/env python3
"""23c_make_arm.py -- write one params.json per history for a clustering arm of the LD-decay test.

Each copy is the derived params.json with r_rate replaced; nothing else changes.
Usage: 23c_make_arm.py <params.json> <arm_dir> <model>=<r_rate> [<model>=<r_rate> ...]
       a bare <r_rate> applies to all seven admissible histories
"""
import json, os, sys

MODELS = ("M0", "M0f", "M1f", "M2", "M2f", "M3", "M3f")
src, arm = sys.argv[1], sys.argv[2]
par = json.load(open(src))
spec = {}
for a in sys.argv[3:]:
    if "=" in a:
        m, r = a.split("=")
        spec[m] = float(r)
    else:
        spec.update({m: float(a) for m in MODELS})
for m, r in spec.items():
    d = os.path.join(arm, m)
    os.makedirs(d, exist_ok=True)
    p = dict(par, r_rate=r)
    json.dump(p, open(os.path.join(d, "params.json"), "w"), indent=1)
    print(f"{d}/params.json r_rate={r:g}")
