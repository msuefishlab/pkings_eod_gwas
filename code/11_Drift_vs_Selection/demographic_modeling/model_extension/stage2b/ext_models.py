#!/usr/bin/env python3
"""ext_models.py -- model extension, Stage 2b: the three alternative model families (CB, CBT, CBTM; each unfolded and
folded, the six alternatives fitted next to the drift-test models).

All are dadi two-population models with BP (Bambomo) = population 1, piecewise-constant sizes and instantaneous
changes, mirroring the drift-test conventions of 07_fit_models.py (sizes relative to N_anc, times in 2*N_anc generations,
dadi m12 = migration into pop 1 (BP) from pop 2 (TP) = "m_BP_from_TP", split by phi_1D_to_2D, default grid, linear
extrapolation over pts = [n+10, n+20, n+30], the misidentification wrapper added by 07's model_setup).

  CB    nu_BP, nu_TP, T_pre, T_c, r_BP                       split at T_s = T_pre + T_c before present; TP constant
                                                              nu_TP; BP nu_BP, changing to r_BP*nu_BP at T_c; no migration
  CBT   + r_TP                                               TP also changes nu_TP -> r_TP*nu_TP at the same T_c
  CBTM  + m_BP_from_TP, m_TP_from_BP, f_mig                   CBT plus asymmetric migration from the split until
                                                              T_m = f_mig*T_s before present (0 <= T_m <= T_s); none after

Reparameterization: the natural parameters are nu_B1, nu_B2, nu_T, T_s, T_c (+ nu_T2; + m_BT, m_TB, T_m).
Here nu_B2 = r_BP*nu_BP, nu_T2 = r_TP*nu_TP, T_s = T_pre + T_c, T_m = f_mig*T_s, so that the contraction-ratio bounds
[1/50, 1], T_c < T_s and 0 <= T_m <= T_s are all box bounds. Nesting values: r_BP = 1 (CB -> M0), r_TP = 1 (CBT -> CB),
m = 0 (CBTM -> CBT; also f_mig = 1), r_BP = r_TP = 1 and f_mig = 0 (CBTM -> M2, and -> M1 with equal m).

Consecutive epochs with identical (nu1, nu2, m12, m21) are merged into one dadi integration, so a model at its nesting
values performs exactly the integration of the nested model (dadi's time stepping depends on how an epoch is split).

T_c window: with N_anc = theta/(4 mu L) and T_c,gen = 2 N_anc T_c, the rule
10^2 <= T_c,gen(mu = 4e-9) and T_c,gen(mu = 1e-8) <= 10^4 is  x = theta*T_c/(2L) in [4e-7, 1e-4]
(wide-bound sensitivity: 10 and 3e4 generations -> x in [4e-8, 3e-4]). theta is not a free dadi parameter, so
the window is applied as a box on T_c computed from a reference theta (tc_box), and checked exactly post-fit.
"""
import math

import numpy as np
from dadi import Integration, Numerics, PhiManip, Spectrum

MERGE_EPOCHS = True

PNAMES = {
    "CB": ["nu_BP", "nu_TP", "T_pre", "T_c", "r_BP"],
    "CBT": ["nu_BP", "nu_TP", "T_pre", "T_c", "r_BP", "r_TP"],
    "CBTM": ["nu_BP", "nu_TP", "T_pre", "T_c", "r_BP", "r_TP", "m_BP_from_TP", "m_TP_from_BP", "f_mig"],
}
X_WINDOW = {"primary": (100 * 4e-9, 1e4 * 1e-8), "wide": (10 * 4e-9, 3e4 * 1e-8)}   # mu * T_c,gen
R_LO = {"primary": 1 / 50, "wide": 1 / 200}
MU = (4e-9, 1e-8)


def tc_box(theta, L, mode="primary"):
    """Box on T_c (dadi units) that maps the generation window onto theta: T_c = 2 L x / theta."""
    x_lo, x_hi = X_WINDOW[mode]
    return 2 * L * x_lo / theta, 2 * L * x_hi / theta


def tc_x(theta, T_c, L):
    """x = mu * T_c,gen = theta * T_c / (2 L); the window is on x."""
    return theta * T_c / (2 * L)


def epochs(nuB, nuT, T_pre, T_c, rB, rT, m12, m21, f_mig):
    """Forward-time epochs [duration, nu1, nu2, m12, m21] from the split (t = 0) to the present (t = T_s)."""
    T_pre, T_c = float(T_pre), float(T_c)
    T_s = T_pre + T_c
    f = min(max(float(f_mig), 0.0), 1.0)
    t_c = min(max(T_pre, 0.0), T_s)            # forward time of the size change
    t_m = min(max(T_s - f * T_s, 0.0), T_s)    # forward time at which migration stops
    cuts = sorted({0.0, t_c, t_m, T_s})
    out = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        if b - a <= 0:
            continue
        mid = 0.5 * (a + b)
        n1, n2 = (float(nuB), float(nuT)) if mid < t_c else (float(rB) * float(nuB), float(rT) * float(nuT))
        mm = (float(m12), float(m21)) if mid < t_m else (0.0, 0.0)
        ep = [b - a, n1, n2, mm[0], mm[1]]
        if MERGE_EPOCHS and out and out[-1][1:] == ep[1:]:
            out[-1][0] += ep[0]
        else:
            out.append(ep)
    return out


def _integrate(eps, ns, pts):
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    for T, n1, n2, a, b in eps:
        phi = Integration.two_pops(phi, xx, T, n1, n2, m12=a, m21=b)
    return Spectrum.from_phi(phi, ns, (xx, xx))


def CB(params, ns, pts):
    nuB, nuT, T_pre, T_c, rB = params
    return _integrate(epochs(nuB, nuT, T_pre, T_c, rB, 1.0, 0.0, 0.0, 1.0), ns, pts)


def CBT(params, ns, pts):
    nuB, nuT, T_pre, T_c, rB, rT = params
    return _integrate(epochs(nuB, nuT, T_pre, T_c, rB, rT, 0.0, 0.0, 1.0), ns, pts)


def CBTM(params, ns, pts):
    nuB, nuT, T_pre, T_c, rB, rT, m12, m21, f_mig = params
    return _integrate(epochs(nuB, nuT, T_pre, T_c, rB, rT, m12, m21, f_mig), ns, pts)


FUNCS = {"CB": CB, "CBT": CBT, "CBTM": CBTM}
EXT = set(FUNCS) | {m + "f" for m in FUNCS}
# nested sub-models (same folding) for the nested-start guard and nested_ok: the chain M0 c CB c CBT c CBTM,
# plus M1, M2 c CBTM (declared extension: CBTM at r = 1, f_mig = 0 is M2)
SUBMODELS = {"CB": ["M0"], "CBT": ["M0", "CB"], "CBTM": ["M0", "CB", "CBT", "M1", "M2"]}


def register(fitmod, mode="primary", theta_ref=None, L=None):
    """Add the three models to an (unchanged) 07_fit_models module's MODELS / BOUNDS / START dicts."""
    for k, f in FUNCS.items():
        fitmod.MODELS[k] = (f, list(PNAMES[k]), None)
    r_lo = R_LO[mode]
    # T_pre carries the drift test's T bounds; migration bounds already in BOUNDS (drift-test M2/M3: 0-30)
    fitmod.BOUNDS.update({"T_pre": (1e-5, 10), "r_BP": (r_lo, 1.0), "r_TP": (r_lo, 1.0), "f_mig": (0.0, 1.0)})
    if theta_ref is not None:
        lo, hi = tc_box(theta_ref, L, mode)
    else:                      # collection / Godambe only: bounds are not used there
        lo, hi = 1e-5, 10.0
    fitmod.BOUNDS["T_c"] = (lo, hi)
    # starts: drift-test values where the parameter exists there (nu 1, T 0.2, m 1); otherwise the centre of the allowed
    # range (geometric for ratio and T_c, arithmetic for the fraction f_mig)
    fitmod.START.update({"T_pre": 0.2, "T_c": math.sqrt(lo * hi), "r_BP": math.sqrt(r_lo), "r_TP": math.sqrt(r_lo),
                         "f_mig": 0.5})
    return lo, hi


def derived(model, p):
    """dadi-unit derived quantities from a parameter dict (base model name, no 'f')."""
    b = model.rstrip("f")
    T_s = p["T_pre"] + p["T_c"]
    d = {"nu_BP2": p["r_BP"] * p["nu_BP"], "nu_TP2": p.get("r_TP", 1.0) * p["nu_TP"], "T_s": T_s}
    if b == "CBTM":
        d["T_m"] = p["f_mig"] * T_s
    return d


def physical(model, p, theta, L, mu):
    """Physical units (simulation/01_fitted_history.py conversion): N_anc = theta/(4 mu L), N = nu N_anc,
    T_gen = 2 N_anc T, m/gen = m/(2 N_anc)."""
    b = model.rstrip("f")
    Na = theta / (4 * mu * L)
    d = derived(model, p)
    out = {"N_anc": Na, "N_BP_1": p["nu_BP"] * Na, "N_BP_2": d["nu_BP2"] * Na, "N_TP_1": p["nu_TP"] * Na,
           "N_TP_2": d["nu_TP2"] * Na, "T_s_gen": 2 * Na * d["T_s"], "T_c_gen": 2 * Na * p["T_c"]}
    if b == "CBTM":
        out["T_m_gen"] = 2 * Na * d["T_m"]
        out["m_BP_from_TP_pergen"] = p["m_BP_from_TP"] / (2 * Na)
        out["m_TP_from_BP_pergen"] = p["m_TP_from_BP"] / (2 * Na)
    return out
