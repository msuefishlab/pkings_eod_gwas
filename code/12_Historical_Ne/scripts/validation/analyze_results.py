#!/usr/bin/env python3
"""
Parse smc++ model.final.json for validation summary and figures.
T_b is estimated as the transition time from "bottleneck-like" Ne to
"ancient-like" Ne, using the geometric-mean of N_recent and N_ancient
as the threshold.
"""
import os, sys, json, csv, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

MU     = 1e-8
N_ANCT = 47250
N_RECT = 2500
# geometric-mean threshold between recent and ancient Ne
N_THRESH = np.sqrt(N_ANCT * N_RECT)   # ≈ 10,869

SCENARIOS = [
    ("null",      None,   N_ANCT, "Constant Ne (null)"),
    ("geneflow",  None,   N_ANCT, "Constant Ne + gene flow (null)"),
    ("Tb500",     500,    N_RECT, "Bottleneck T_b=500 gen"),
    ("Tb1000",    1000,   N_RECT, "Bottleneck T_b=1000 gen"),
    ("Tb2000",    2000,   N_RECT, "Bottleneck T_b=2000 gen"),
    ("Tb5000",    5000,   N_RECT, "Bottleneck T_b=5000 gen"),
    ("Tb10000",   10000,  N_RECT, "Bottleneck T_b=10000 gen"),
]

def smcpp_json_to_steps(json_file, mu=1e-8):
    """
    Parse smc++ model.final.json.
    Returns (t_knots_gen, Ne_at_knot, N0) for a piecewise-constant model.
    t_knots_gen[i] is the BEGINNING of interval i (in generations).
    Ne_at_knot[i] is Ne during interval [t_knots_gen[i], t_knots_gen[i+1]).
    """
    with open(json_file) as f:
        d = json.load(f)
    theta  = d["theta"]
    model  = d["model"]
    N0     = theta / (2.0 * mu)
    knots  = np.array(model["knots"])
    y      = np.array(model["y"])
    t_gen  = knots * 2.0 * N0
    Ne     = N0 * np.exp(y)
    return t_gen, Ne, N0

def infer_Tb_and_depth(t_gen, Ne, thresh=N_THRESH):
    """
    For a piecewise-constant Ne(t):
    T_b = earliest time where Ne FIRST exceeds thresh (transition onset).
    Ne_min = minimum Ne across all knots.
    """
    Ne_min = float(np.min(Ne))
    # Find the first knot where Ne jumps above the threshold
    for i in range(len(Ne)):
        if Ne[i] > thresh:
            Tb_inf = float(t_gen[i])
            return Tb_inf, Ne_min
    return None, Ne_min  # never crosses threshold

def true_ne_curve(scenario, true_T_b, true_N_recent, t_max=200_000):
    t = np.logspace(np.log10(50), np.log10(t_max), 400)
    if scenario in ("null", "geneflow"):
        ne = np.full_like(t, N_ANCT, dtype=float)
    else:
        ne = np.where(t < true_T_b, float(true_N_recent), float(N_ANCT))
    return t, ne

def piecewise_ne_curve(t_gen, Ne, t_max=200_000):
    """Return dense (t, Ne) arrays from piecewise-constant steps."""
    t = np.logspace(np.log10(max(50, t_gen[0] * 0.3)),
                    np.log10(t_max), 400)
    ne = np.empty_like(t)
    for j, tj in enumerate(t):
        idx = np.searchsorted(t_gen, tj, side="right") - 1
        idx = max(0, min(idx, len(Ne) - 1))
        ne[j] = Ne[idx]
    return t, ne

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir",   required=True)
    ap.add_argument("--outdir",    required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    summary_rows = []
    summary_path = os.path.join(args.outdir, "validation_summary.tsv")

    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    axes = axes.flatten()
    axes[7].set_visible(False)

    for i, (sc, true_T_b, true_N_recent, label) in enumerate(SCENARIOS):
        ax = axes[i]
        json_f = os.path.join(args.workdir, "estimates", sc, "model.final.json")

        if not os.path.exists(json_f):
            ax.set_title(f"{label}\n[MISSING]", fontsize=8)
            ax.text(0.5, 0.5, "No output", ha="center", va="center",
                    transform=ax.transAxes, color="red")
            summary_rows.append({
                "scenario": sc, "true_T_b": "NA" if true_T_b is None else str(true_T_b),
                "true_N_recent": str(true_N_recent), "inferred_T_b": "NA",
                "inferred_N_recent": "NA", "recovered": "NA", "false_positive": "NA",
                "notes": "MISSING"})
            continue

        t_gen, Ne_step, N0 = smcpp_json_to_steps(json_f, mu=MU)
        inf_Tb, inf_Ne_min = infer_Tb_and_depth(t_gen, Ne_step)
        t_inf, ne_inf = piecewise_ne_curve(t_gen, Ne_step)
        t_true, ne_true = true_ne_curve(sc, true_T_b, true_N_recent)

        print(f"  {sc}: N0={N0:.0f}", flush=True)
        print(f"    knots(gen): {list(t_gen.astype(int))}", flush=True)
        print(f"    Ne at knots: {[f'{n:.0f}' for n in Ne_step]}", flush=True)
        tb_s  = f"{inf_Tb:.0f}"  if inf_Tb  is not None else "NA"
        ne_s  = f"{inf_Ne_min:.0f}" if inf_Ne_min is not None else "NA"
        print(f"    inferred T_b={tb_s}, Ne_min={ne_s}", flush=True)

        # Evaluate recovery
        notes = ""
        if true_T_b is not None:
            if inf_Tb is not None:
                t_ratio = inf_Tb / true_T_b
                t_ok = 0.5 <= t_ratio <= 2.0
            else:
                t_ok = False
                t_ratio = None
            ne_ratio = inf_Ne_min / true_N_recent if inf_Ne_min else None
            # depth: allow factor of 3 (systematic downward bias with r/mu=1 data)
            ne_ok = (ne_ratio is not None and 1/3 <= ne_ratio <= 3.0)

            if t_ok and ne_ok:
                recovered = "YES"
            elif t_ok:
                recovered = "PARTIAL (depth off)"
                notes = f"depth_ratio={ne_ratio:.2f}"
            elif ne_ok:
                recovered = "PARTIAL (timing off)"
                notes = f"timing_ratio={t_ratio:.2f}" if t_ratio else "Tb not detected"
            else:
                recovered = "NO"
                notes = (f"timing_ratio={t_ratio:.2f}, " if t_ratio else "") + \
                        (f"depth_ratio={ne_ratio:.2f}" if ne_ratio else "")
            fp = "NA"
        else:
            recovered = "NA"
            # false positive: inferred minimum Ne < N_ancient/3 at recent times (< 50k gen)
            mask_recent = t_inf < 50_000
            ne_min_recent = float(np.min(ne_inf[mask_recent])) if np.any(mask_recent) else None
            fp = ("YES" if (ne_min_recent is not None and ne_min_recent < N_ANCT / 3)
                  else "NO")
            notes = f"Ne_min_recent={ne_min_recent:.0f}" if ne_min_recent else ""

        summary_rows.append({
            "scenario":          sc,
            "true_T_b":          "NA" if true_T_b is None else str(true_T_b),
            "true_N_recent":     str(true_N_recent),
            "inferred_T_b":      "NA" if inf_Tb is None else f"{inf_Tb:.0f}",
            "inferred_N_recent": "NA" if inf_Ne_min is None else f"{inf_Ne_min:.0f}",
            "recovered":         recovered,
            "false_positive":    fp,
            "notes":             notes,
        })

        # plot
        ax.plot(t_true, ne_true, "k--", lw=1.8, alpha=0.7, label="True", zorder=5)
        # step function (the actual piecewise model)
        t_step = np.concatenate([[t_gen[0] * 0.3], np.repeat(t_gen[1:], 2),
                                  [t_gen[-1] * 3]])
        ne_step_rep = np.repeat(Ne_step, 2)
        ax.plot(t_step, ne_step_rep, color="steelblue", lw=1.5,
                label="smc++ (piecewise)", zorder=4)
        if true_T_b:
            ax.axvline(true_T_b, color="gray", ls=":", lw=1.2, alpha=0.7,
                       label=f"True T_b={true_T_b:,}")
        if inf_Tb:
            ax.axvline(inf_Tb, color="tomato", ls="--", lw=1.2, alpha=0.7,
                       label=f"Inferred T_b={inf_Tb:.0f}")
        ax.axhline(N_THRESH, color="gray", ls="-.", lw=0.8, alpha=0.5)

        ax.set_xscale("log"); ax.set_yscale("log")
        status = f"recovered={recovered}  FP={fp}" if sc in ("null", "geneflow") else \
                 f"recovered={recovered}"
        ax.set_title(f"{label}\n{status}", fontsize=7)
        ax.set_xlabel("Time (gen)", fontsize=7)
        ax.set_ylabel("Ne", fontsize=7)
        ax.legend(fontsize=6, loc="upper right")
        ax.tick_params(labelsize=7)
        ax.set_xlim(50, 300_000)
        ax.set_ylim(300, 2_000_000)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: f"{int(x):,}"))

    fig.suptitle(
        "smc++ recovery validation — Pkingsleyae parameters\n"
        "msprime: 20 diploids, 25 chr (~934 Mb), μ=r=1×10⁻⁸, 20% masked | "
        "N_ancient=47,250, N_recent=2,500",
        fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    fig_path = os.path.join(args.outdir, "validation_figure.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nFigure saved: {fig_path}", flush=True)

    # write TSV
    fields = ["scenario", "true_T_b", "true_N_recent",
              "inferred_T_b", "inferred_N_recent", "recovered",
              "false_positive", "notes"]
    with open(summary_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)
    print(f"Summary saved: {summary_path}", flush=True)

    print("\n" + "="*70)
    print("VALIDATION SUMMARY")
    print("="*70)
    print(f"{'scenario':<12} {'true_T_b':>8} {'inf_T_b':>9} "
          f"{'true_Ne':>8} {'inf_Ne':>8}  recovered    FP")
    print("-" * 70)
    for row in summary_rows:
        print(f"{row['scenario']:<12} {row['true_T_b']:>8} "
              f"{row['inferred_T_b']:>9} {row['true_N_recent']:>8} "
              f"{row['inferred_N_recent']:>8}  "
              f"{row['recovered']:<20} {row['false_positive']}")

if __name__ == "__main__":
    main()
