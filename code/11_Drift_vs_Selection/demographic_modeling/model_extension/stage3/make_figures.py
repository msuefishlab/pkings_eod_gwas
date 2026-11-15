"""make_figures.py -- Stage 3 figure functions (imported after analysis_step5.py; matplotlib style set by the caller).
fig_ld(...)    observed vs simulated binned r^2 (BAM): calibration and held-out panels, primary and callable-only rows
fig_diag(...)  differentiation / supporting diagnostics: replicate distributions per model vs the observed value
"""
import numpy as np
import matplotlib.pyplot as plt

MID = (np.arange(50) + 0.5) * 2


def fig_ld(obs, sims, models, colors, fn, title_rows=("All observed SNVs (primary, gates P1)",
                                                         "Callable SNVs on both sides (sensitivity)")):
    """obs[(row, set)] = dict(r2, lo, hi); sims[(row, set, model)] = curve"""
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.4), sharex=True, sharey=True)
    for i, row in enumerate(("primary", "callable")):
        for j, s in enumerate(("calibration", "heldout")):
            ax = axes[i, j]
            o = obs[(row, s)]
            ax.fill_between(MID, o["lo"], o["hi"], color="0.82", lw=0, label="observed 95% band")
            ax.plot(MID, o["r2"], color="k", lw=1.4, label="observed")
            for m in models:
                c = sims.get((row, s, m))
                if c is not None:
                    ax.plot(MID, c, color=colors[m], lw=1.0, label=m)
            ax.set_title(f"{title_rows[i]}\n{'calibration chromosomes (fit)' if s == 'calibration' else 'held-out chromosomes (scored)'}",
                         loc="left")
            if i == 1:
                ax.set_xlabel("distance (kb, callable-layout units)")
            if j == 0:
                ax.set_ylabel("mean r² (BAM)")
            ax.margins(x=0.02)
    axes[0, 1].legend(frameon=False, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(fn, dpi=300)
    return fig


def fig_diag(stats, obsv, dists, models, colors, fn):
    """stats: list of (key, label); obsv[key] observed value; dists[(key, model)] = array of replicate values"""
    k = len(stats)
    ncol = 4
    nrow = int(np.ceil(k / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2, 1.9 * nrow + 0.4))
    axes = np.atleast_2d(axes)
    for a, (key, lab) in enumerate(stats):
        ax = axes[a // ncol, a % ncol]
        for i, m in enumerate(models):
            v = dists.get((key, m))
            if v is None or len(v) == 0:
                continue
            lo, hi = np.percentile(v, [2.5, 97.5])
            ax.plot([lo, hi], [i, i], color=colors[m], lw=2.2, solid_capstyle="butt")
            ax.plot([np.median(v)], [i], "o", color=colors[m], ms=3)
        ax.axvline(obsv[key], color="k", lw=1.0)
        ax.set_yticks(range(len(models)))
        ax.set_yticklabels(models if a % ncol == 0 else [])
        ax.set_title(lab, loc="left")
        ax.margins(x=0.08, y=0.12)
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, 4))
    for a in range(k, nrow * ncol):
        axes[a // ncol, a % ncol].axis("off")
    fig.tight_layout()
    fig.savefig(fn, dpi=300)
    return fig
