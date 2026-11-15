#!/usr/bin/env python3
"""Chromosome-block bootstrap dirs for the v2 tree: replicate K resamples the 25
chromosomes with replacement (seed K*1000 + pop_index*10 + mu_index) and links all
5 DIs of each drawn chromosome. Usage: make_bootstraps.py <v2_tree> <manifest>.
Writes tasks (TAB) to argv[2]."""
import os, random, sys
V2 = sys.argv[1]; out_manifest = sys.argv[2]
POPS = ["BAM", "BIK", "BAVA", "BIR", "APA_2009", "BEN", "DOV", "DOG"]
MU_ARMS = [("1e-8", "1e8"), ("4e-9", "4e9")]
CHR = [f"chr{i}" for i in range(1, 26)]; DIS = [f"di{j}" for j in range(1, 6)]
tasks, missing = [], 0
for K in range(1, 11):
    for pi, pop in enumerate(POPS):
        if not all(os.path.exists(f"{V2}/smc_input/{pop}/{c}_{d}.smc.gz") for c in CHR for d in DIS):
            print(f"skip {pop}: marginal inputs incomplete"); continue
        for mi, (mu, lab) in enumerate(MU_ARMS):
            rng = random.Random(K * 1000 + pi * 10 + mi)
            draw = rng.choices(CHR, k=25)
            bd = f"{V2}/bootstraps/boot_{K}/{pop}_mu{lab}"; os.makedirs(bd, exist_ok=True)
            for i, oc in enumerate(draw, start=1):
                for d in DIS:
                    dst = f"{bd}/chr{i}_{d}.smc.gz"
                    if os.path.lexists(dst): os.remove(dst)
                    os.symlink(f"{V2}/smc_input/{pop}/{oc}_{d}.smc.gz", dst)
            tasks.append(f"boot\t{K}\t{pop}\t{mu}\t{lab}")
with open(out_manifest, "w") as f:
    f.write("".join(t + "\n" for t in tasks))
print(f"{len(tasks)} bootstrap tasks -> {out_manifest}")
