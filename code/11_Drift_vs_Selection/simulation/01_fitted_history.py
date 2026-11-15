#!/usr/bin/env python3
"""01_fitted_history.py -- build the simulation's params.json from the drift-test outputs
(demographic_modeling/fixed_differences/).

The fitted A1w histories are READ FROM THE DRIFT TEST, not from a scratch copy: everything below
comes out of output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences{,_agree}/.

  tail/models.tsv        per model: theta, T (2N_anc units), nu_BP, nu_TP, m_*, s_BP, N_anc per mu,
                         and the adequate / admissible flags
  checks/spectra_checks.tsv   L_mask (check "total") and L_mask_alt (check "total_alt")
  sfs/<pair>.block_scan.tsv   L_scan / L_scan_alt, summed over blocks
  tail/index_units.tsv        the observed Peak 4 / Peak 5 delta, from the integer allele counts

Unit conversion, exactly as METHODS.md section 2 states it:
  N_anc  = theta / (4 mu L_mask)            (folded models use L_mask_alt)
  N_i    = nu_i * N_anc                     (M3/M3f: start = s_BP * N_anc / (1 - s_BP) * N_anc)
  T_gen  = T * 2 N_anc
  m_pergen = m_dadi / (2 N_anc)             (dadi naming: m_BP_from_TP is forwards-time)

Two things in params.json are NOT fitted and are declared assumptions, carried here as constants
with their METHODS citation: the recombination rate r = 1e-8 /bp/gen (METHODS.md section 2 --
there is no recombination map for this species) and the 2-Mb simulation chunk length.

Usage:
  01_fitted_history.py --agree-dir <fd_out> --out <params.json> [--pair A1w] [--select admissible]
  01_fitted_history.py --agree-dir <fd_out> --check params.as_run.agree.json
"""
import argparse, csv, json, math, os, sys

MU = (4e-9, 1e-8)          # 09_tail.py line 61
R_RATE = 1e-8              # METHODS.md section 2: chosen, not inherited. No map exists.
CHUNK_LEN = 2_000_000
SIGFIG = 6
# Keys an as-run params.json may carry that are NOT derivable from the drift-test outputs and are
# not read by any worker; absent from the derived file by design.
NON_DERIVED = ("GENOME",)


def sig(x, n=SIGFIG):
    if x == 0:
        return 0.0
    return float("%.*g" % (n, x))


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_params(s):
    out = {}
    for kv in s.split(";"):
        if not kv:
            continue
        k, v = kv.split("=")
        out[k] = float(v)
    return out


def lengths(fd_out, pair):
    """L_mask / L_mask_alt from the spectra gate, L_scan / L_scan_alt from the block scan."""
    L = {}
    for r in read_tsv(os.path.join(fd_out, "checks", "spectra_checks.tsv")):
        if r["pair"] != pair:
            continue
        if r["check"] == "total":
            L["L_mask"] = float(r["value"])
        elif r["check"] == "total_alt":
            L["L_mask_alt"] = float(r["value"])
    scan = scan_alt = 0.0
    for r in read_tsv(os.path.join(fd_out, "sfs", "%s.block_scan.tsv" % pair)):
        scan += float(r["L_scan"])
        scan_alt += float(r.get("L_scan_alt") or 0.0)
    L["L_scan"], L["L_scan_alt"] = scan, scan_alt
    if "L_mask" not in L:
        sys.exit("ERROR: L_mask (check \"total\") not in checks/spectra_checks.tsv for pair %s" % pair)
    return L   # L_mask_alt / L_scan_alt may be absent when the run has no folded models


def deltas(fd_out, pair):
    """Observed Peak 4 / Peak 5 delta, recomputed from the integer allele counts."""
    out = {}
    for r in read_tsv(os.path.join(fd_out, "tail", "index_units.tsv")):
        if r["pair_id"] != pair or r["peak"] not in ("4", "5"):
            continue
        out["delta_peak%s" % r["peak"]] = (int(r["d_bp"]) / int(r["AN_bp"])
                                          - int(r["d_tp"]) / int(r["AN_tp"]))
    for k in ("delta_peak4", "delta_peak5"):
        if k not in out:
            sys.exit("ERROR: %s not derivable from tail/index_units.tsv" % k)
    return out


def build(fd_out, pair, select):
    L = lengths(fd_out, pair)
    rows = [r for r in read_tsv(os.path.join(fd_out, "tail", "models.tsv")) if r["pair"] == pair]
    if not rows:
        sys.exit("ERROR: pair %s not in tail/models.tsv" % pair)
    kept = [r for r in rows if r[select] == "True" and r["converged"] == "True"]
    if not kept:
        sys.exit("ERROR: no model for %s passes converged + %s" % (pair, select))
    par = {}
    for r in kept:
        model = r["model"]
        folded = r["analysis"] == "folded"
        if folded and "L_mask_alt" not in L:
            sys.exit("ERROR: model %s is folded but checks/spectra_checks.tsv has no \"total_alt\" "
                     "row for pair %s" % (model, pair))
        l_mask = L["L_mask_alt"] if folded else L["L_mask"]
        l_scan = L["L_scan_alt"] if folded else L["L_scan"]
        p = parse_params(r["params"])
        theta = float(r["theta"])
        for mu in MU:
            n_anc = theta / (4.0 * mu * l_mask)
            bp_end = p["nu_BP"] * n_anc
            tp_end = p["nu_TP"] * n_anc
            if "s_BP" in p:                       # M3 / M3f: exponential growth from the split
                bp_start = p["s_BP"] * n_anc
                tp_start = (1.0 - p["s_BP"]) * n_anc
            else:
                bp_start, tp_start = bp_end, tp_end
            m_bp = p.get("m_BP_from_TP", p.get("m", 0.0))
            m_tp = p.get("m_TP_from_BP", p.get("m", 0.0))
            par["%s|%g" % (model, mu)] = {
                "N_anc": sig(n_anc),
                "N_BP_start": sig(bp_start), "N_BP_end": sig(bp_end),
                "N_TP_start": sig(tp_start), "N_TP_end": sig(tp_end),
                "T_gen": sig(p["T"] * 2.0 * n_anc),
                "m_BP_from_TP_pergen": sig(m_bp / (2.0 * n_anc)),
                "m_TP_from_BP_pergen": sig(m_tp / (2.0 * n_anc)),
                "L_mask": sig(l_mask), "L_scan": sig(l_scan),
            }
    par["r_rate"] = R_RATE
    par["chunk_len"] = CHUNK_LEN
    par.update(deltas(fd_out, pair))
    par["L_mask"] = sig(L["L_mask"])
    par["L_scan"] = sig(L["L_scan"])
    return par, sorted(r["model"] for r in kept)


def compare(new, ref_path, tol):
    ref = json.load(open(ref_path))
    only_ref = sorted(set(ref) - set(new) - set(NON_DERIVED))
    only_new = sorted(set(new) - set(ref))
    worst, worst_key = 0.0, None
    for k in sorted(set(ref) & set(new)):
        a, b = ref[k], new[k]
        items = [(k, a, b)] if not isinstance(a, dict) else [
            ("%s.%s" % (k, kk), a[kk], b[kk]) for kk in a if kk in b]
        for name, x, y in items:
            if x == 0 and y == 0:
                continue
            rel = abs(y - x) / max(abs(x), 1e-300)
            if rel > worst:
                worst, worst_key = rel, name
    skipped = sorted(set(ref) & set(NON_DERIVED))
    print("reference        : %s" % ref_path)
    print("not derivable    : %s" % (", ".join(skipped) or "-"))
    print("keys only in ref : %s" % (", ".join(only_ref) or "-"))
    print("keys only in new : %s" % (", ".join(only_new) or "-"))
    print("max rel. diff    : %.3g  (%s)" % (worst, worst_key))
    print("tolerance        : %.3g" % tol)
    if only_ref:
        sys.exit("FAIL: the derivation does not cover every as-run model")
    if worst > tol:
        sys.exit("FAIL: derived history differs from the as-run history by more than the tolerance")
    print("PASS")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--agree-dir", required=True,
                    help="output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree")
    ap.add_argument("--pair", default="A1w")
    ap.add_argument("--select", default="admissible", choices=["admissible", "adequate"],
                    help="admissible = the agreement run's 7 models; adequate = the first-pass A1w "
                         "run, where no model is admissible (see README.md)")
    ap.add_argument("--out")
    ap.add_argument("--check", help="compare against an as-run params.json instead of writing")
    ap.add_argument("--tol", type=float, default=1e-4)
    a = ap.parse_args()
    par, models = build(a.agree_dir, a.pair, a.select)
    print("pair %s, %s models: %s" % (a.pair, a.select, ", ".join(models)), file=sys.stderr)
    if a.check:
        compare(par, a.check, a.tol)
        return
    if not a.out:
        sys.exit("ERROR: --out or --check required")
    with open(a.out, "w") as fh:
        json.dump(par, fh, indent=1)
    print(a.out)


if __name__ == "__main__":
    main()
