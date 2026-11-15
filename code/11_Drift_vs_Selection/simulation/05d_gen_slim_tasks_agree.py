#!/usr/bin/env python3
import sys, random
Q = 10
MODELS = {'M0': (12184, 30588, 0.0, 0.0, 4348), 'M0f': (15011, 40156, 0.0, 0.0, 5159), 'M1f': (18448, 45077, 2.96885e-05, 2.96885e-05, 248448), 'M2': (17499, 32770, 2.49228e-05, 4.13148e-05, 170103), 'M2f': (21186, 41221, 2.02159e-05, 3.91774e-05, 203358), 'M3': (17672, 31919, 2.43284e-05, 4.23664e-05, 428591), 'M3f': (21288, 36220, 2.4011e-05, 3.91257e-05, 121067)}
H_VALS = [0.0, 0.5, 1.0]
S_VALS = [0.0002, 0.0005, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.01, 0.02, 0.05, 0.1]
ONSETS = [25, 50, 100, 250, 500, 1000, 2500, 5000]
F0 = {"denovo": 0.0, "st001": 0.01, "st005": 0.05}
REPS = {"denovo": 100, "st001": 50, "st005": 50}
rng = random.Random(20261022)
out = open(sys.argv[1], "w"); n = 0
def emit(M, h, s, onset, f0, seltpf, reps):
    global n
    nbp, ntp, mbp, mtp, tg = MODELS[M]
    if onset > tg:
        return
    nbpq, ntpq = max(100, round(nbp / Q)), max(100, round(ntp / Q))
    sq = s * Q; tend = max(2, round(onset / Q)) + 1
    seltpq = -sq if seltpf else 0.0
    for rep in range(1, reps + 1):
        out.write(f"{M} {Q} {onset} {seltpf} {nbpq} {ntpq} {mbp*Q:g} {mtp*Q:g} {sq:g} {h:g} {f0:g} "
                  f"{tend} {seltpq:g} {rep} {rng.randrange(1,2_000_000_000)} {s:g}\n")
        n += 1
for M in MODELS:
    for f0n, f0 in F0.items():
        for onset in ONSETS:
            emit(M, 0.5, 0.0, onset, f0, 0, 500)          # neutral controls
    for h in H_VALS:
        for s in S_VALS:
            for onset in ONSETS:
                for f0n, f0 in F0.items():
                    for seltpf in (0, 1):
                        emit(M, h, s, onset, f0, seltpf, REPS[f0n])
out.close(); print(n)
