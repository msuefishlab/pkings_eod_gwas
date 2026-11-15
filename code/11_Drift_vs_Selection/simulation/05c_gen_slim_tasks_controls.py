#!/usr/bin/env python3
import sys, random
Q = 10
MODELS = {"M0": (14205, 38898, 0.0, 0.0, 5104), "M1": (17734, 43126, 3.0e-5, 3.0e-5, 283387),
          "M2": (20075, 40227, 2.2e-5, 3.7e-5, 69813), "M3": (20097, 39419, 2.1e-5, 3.7e-5, 341605)}
rng = random.Random(20261021)
out = open(sys.argv[1], "w"); n = 0
def emit(M, h, s, onset, f0, seltpf, reps):
    global n
    nbp, ntp, mbp, mtp, tg = MODELS[M]
    if onset > tg: return
    nbpq, ntpq = max(100, round(nbp / Q)), max(100, round(ntp / Q))
    sq = s * Q; tend = max(2, round(onset / Q)) + 1
    seltpq = -sq if seltpf else 0.0
    for rep in range(1, reps + 1):
        out.write(f"{M} {Q} {onset} {seltpf} {nbpq} {ntpq} {mbp*Q:g} {mtp*Q:g} {sq:g} {h:g} {f0:g} "
                  f"{tend} {seltpq:g} {rep+100000} {rng.randrange(1,2_000_000_000)} {s:g}\n")
        n += 1
# (a) neutral controls: s = 0, every model x f0 x onset, 500 replicates
for M in MODELS:
    for f0 in (0.0, 0.01, 0.05):
        for onset in (25, 50, 100, 250, 500, 1000, 2500, 5000):
            emit(M, 0.5, 0.0, onset, f0, 0, 500)
# (b) finer s near the power floor, at the onsets where the main grid found compatible cells
for M in MODELS:
    for h in (0.0, 0.5, 1.0):
        for s in (0.0002, 0.0005, 0.001, 0.0015, 0.002, 0.003):
            for onset in (1000, 2500, 5000):
                for f0 in (0.01, 0.05):
                    for seltpf in (0, 1):
                        emit(M, h, s, onset, f0, seltpf, 100)
out.close(); print(n)
