#!/usr/bin/env python3
import sys, random
Q = 10
MODELS = {  # mu = 1e-8 arm: (N_BP_end, N_TP_end, m_BP_from_TP, m_TP_from_BP, T_gen)
 "M0": (14205, 38898, 0.0,     0.0,     5104),
 "M1": (17734, 43126, 3.0e-5,  3.0e-5,  283387),
 "M2": (20075, 40227, 2.2e-5,  3.7e-5,  69813),
 "M3": (20097, 39419, 2.1e-5,  3.7e-5,  341605),
}
H_VALS = [0.0, 0.5, 1.0]
S_VALS = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]
ONSETS = [25, 50, 100, 250, 500, 1000, 2500, 5000]
F0 = {"denovo": 0.0, "st001": 0.01, "st005": 0.05}
REPS = {"denovo": 500, "st001": 200, "st005": 200}
rng = random.Random(20261020)
out = open(sys.argv[1], "w")
n = 0
for M, (nbp, ntp, mbp, mtp, tg) in MODELS.items():
    nbpq, ntpq = max(100, round(nbp / Q)), max(100, round(ntp / Q))
    mbpq, mtpq = mbp * Q, mtp * Q
    for h in H_VALS:
        for s in S_VALS:
            sq = s * Q
            for onset in ONSETS:
                if onset > tg:
                    continue
                tend = max(2, round(onset / Q)) + 1
                for f0name, f0 in F0.items():
                    for seltpf in (0, 1):
                        seltpq = -sq if seltpf else 0.0
                        for rep in range(1, REPS[f0name] + 1):
                            seed = rng.randrange(1, 2_000_000_000)
                            out.write(f"{M} {Q} {onset} {seltpf} {nbpq} {ntpq} {mbpq:g} {mtpq:g} "
                                      f"{sq:g} {h:g} {f0:g} {tend} {seltpq:g} {rep} {seed} {s:g}\n")
                            n += 1
out.close()
print(n)
