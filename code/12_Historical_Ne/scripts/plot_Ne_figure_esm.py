#!/usr/bin/env python3
"""Fig. S9: SMC++ Ne(t) on the corrected masks (masks_v2), mu = 1e-8 arm only.

Panel A, the four BP populations; panel B, the four TP populations. Lines are point estimates and shaded bands the
2.5-97.5 percentiles of the 10 chromosome bootstraps (columns Ne, Ne_lo_95, Ne_hi_95 of collate.py's table).
Okabe-Ito colours plus line styles as a second cue; 8.0 x 3.3 in, fonts >= 7 pt.

Usage (repository root; matplotlib only):
    python3 code/12_Historical_Ne/scripts/plot_Ne_figure_esm.py \
        output_data/12_Historical_Ne/results/Ne_all_pops_both_mu.tsv \
        output_data/12_Historical_Ne/esm_fig_s9/Fig_S9_Ne
"""
import csv
import sys
from collections import defaultdict

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as ticker  # noqa: E402

COLORS = ["#0072B2", "#D55E00", "#009E73", "#E69F00"]
STYLES = ["-", (0, (5, 1.5)), (0, (1.2, 1.2)), (0, (5, 1.5, 1.2, 1.5))]
PANELS = [
    ("A", "Biphasic populations",
     [("BAM", "Bambomo"), ("BIK", "Bikagala"), ("BAVA", "Bavavela"), ("BIR", "Biroundou")]),
    ("B", "Triphasic populations",
     [("APA_2009", "Apassa (2009)"), ("BEN", "Bengue"), ("DOV", "Douvalou"), ("DOG", "Douengi")]),
]
MU = "1e8"


def fmt(x, _):
    if x >= 1e6:
        return "%gM" % (x / 1e6)
    if x >= 1e3:
        return "%gk" % (x / 1e3)
    return "%g" % x


def main(tsv, out_stem):
    data = defaultdict(list)
    with open(tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["mu"] == MU:
                data[r["pop"]].append(tuple(float(r[k]) for k in ("t_gen", "Ne", "Ne_lo_95", "Ne_hi_95")))
    for k in data:
        data[k].sort()

    mpl.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.linewidth": 0.6,
                         "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.4,
                         "ytick.minor.width": 0.4, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.3), sharey=True)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.15, top=0.91, wspace=0.07)
    for ax, (letter, title, pops) in zip(axes, PANELS):
        for (code, name), c, ls in zip(pops, COLORS, STYLES):
            rows = data[code]
            assert rows, code
            t = [r[0] for r in rows]
            ax.fill_between(t, [r[2] for r in rows], [r[3] for r in rows], color=c, alpha=0.13, lw=0)
            ax.plot(t, [r[1] for r in rows], color=c, ls=ls, lw=1.5, label=name)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(10, 2e6)
        ax.set_ylim(500, 3e5)
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(ticker.FuncFormatter(fmt))
            axis.set_minor_formatter(ticker.NullFormatter())
        ax.grid(True, which="major", color="#d0d0d0", lw=0.4)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.set_xlabel("Generations before present", fontsize=8.5)
        ax.set_title(title, fontsize=9, loc="left", pad=4)
        ax.text(-0.02, 1.02, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", ha="right", va="bottom")
        ax.legend(loc="lower right", fontsize=7.5, frameon=True, framealpha=0.9, edgecolor="none",
                  handlelength=2.6, labelspacing=0.3, borderpad=0.4)
    axes[0].set_ylabel(r"Effective population size, $N_e$", fontsize=8.5)
    for ext in ("pdf", "png"):
        fig.savefig("%s.%s" % (out_stem, ext), dpi=300)
        print("saved %s.%s" % (out_stem, ext))


if __name__ == "__main__":
    main(*sys.argv[1:3])
