#!/usr/bin/env python3
"""04_fit_hap.py -- phased two-state local-ancestry fits for every fish outside a panel.

For one marker panel (03_markers.py): per fish, maximum-likelihood (T, m) on the la_lib grid from the
genome-wide markers (both haplotypes share T and m), 95% intervals from 1,000 chromosome bootstrap resamples,
posterior P(BAM) at the fish's ML on (a) the genome-wide markers, for the fraction and tracts, and (b) all
markers incl. peaks and the chr4 inversion, for the peak test (07). T is fitted at r = 1.4 cM/Mb and reported
at 1.0, 1.4 and 3.0 cM/Mb as T(r) = T(1.4) * 1.4 / r.

fits.tsv also carries the statistics used to match fish to simulations (07_match.py): q_m (projection at the AIMs), F_long (fraction of
analysed length in tracts >= 2 Mb) and log10 T_hat (column T_r1.4 is T_hat).

Writes in <out_dir>/<tag>/ (tag = hap_<panel>): fits.tsv, tracts.tsv, surfaces.npy (chrom x T x m x fish),
post_gw.npz and post_all.npz (P(BAM) per marker x haplotype, with sample order).

Usage (repo root): python3 code/11_Drift_vs_Selection/local_ancestry/04_fit_hap.py <out_dir> <panel> [workers]
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import la_lib as L

PANEL_ROLES = {"primary": ("PB", "PA", "BN"), "full": ("FB", "FA", "BN"), "ben": ("PB", "BN", "FA")}  # emission + ascertainment


def read_samples(out):
    info = {}
    for line in open(f"{out}/samples.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] != "sample":
            info[f[0]] = dict(group=f[2], q=f[3], site=f[4], roles=set(f[5].split(",")))
    return info


def main():
    out, panel = sys.argv[1], sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    tag = f"hap_{panel}"
    od = f"{out}/{tag}"; os.makedirs(od, exist_ok=True)
    info = read_samples(out)
    M = L.load_markers(out, panel)
    hs = M["hap_sample"]
    roles = PANEL_ROLES[panel]
    fish = [s for s in dict.fromkeys(hs) if s in info and not (info[s]["roles"] & set(roles))]  # BAM_175 dropped
    cols = np.concatenate([np.flatnonzero(hs == s) for s in fish])
    assert cols.size == 2 * len(fish)
    obs = M["H"][:, cols]
    gw = M["gw"].astype(bool)
    print(f"{tag}: {len(fish)} fish, {gw.sum()} genome-wide markers, {M['pos'].size} all", flush=True)

    S = L.grid_surfaces("hap", L.per_chrom_split(M["chrom"], M["pos"], obs, M["pBe"], M["pTe"], gw), workers)
    Sf = S[..., 0::2] + S[..., 1::2]                      # sum the two haplotypes of each fish
    np.save(f"{od}/surfaces.npy", Sf.astype(np.float32))
    fits = L.fit_from_surfaces(Sf)
    lam = np.repeat([f["T"] * L.R_FIT for f in fits], 2)
    mm = np.repeat([f["m"] for f in fits], 2)

    excl = L.excluded_regions()
    for label, mask in (("gw", gw), ("all", np.ones(gw.size, bool))):
        post = np.empty((mask.sum(), obs.shape[1]))
        idx = np.flatnonzero(mask)
        for c in range(1, 26):
            s = M["chrom"][idx] == c
            ii = idx[s]
            post[s] = L.hap_posterior(obs[ii], M["pos"][ii], M["pBe"][ii], M["pTe"][ii], lam, mm)
        np.savez_compressed(f"{od}/post_{label}.npz", post=post.astype(np.float32), chrom=M["chrom"][idx],
                            pos=M["pos"][idx], peak=M["peak"][idx], hap_sample=np.repeat(fish, 2))
        if label == "gw":
            post_gw = post
    seg = L.segments(M["chrom"][gw], M["pos"][gw], excl)
    rows, trows = [], []
    for i, s in enumerate(fish):
        f = fits[i]
        sm, tr = L.summarise_fish(post_gw[:, 2 * i:2 * i + 2], M["pos"][gw], M["chrom"][gw], seg, f["T"], f["m"])
        for h, ch, a, b in tr:
            trows.append(f"{s}\t{h}\tchr{ch}\t{a:.0f}\t{b:.0f}\t{b - a:.0f}")
        Tr = {r: (f["T"] * 1.4 / r, f["T_lo"] * 1.4 / r, f["T_hi"] * 1.4 / r) for r in L.R_BRACKET}
        rows.append("\t".join([s, info[s]["group"], info[s]["site"], info[s]["q"],
                               *(f"{x:.1f}" for r in L.R_BRACKET for x in Tr[r]),
                               f"{f['m']:.3f}", f"{f['m_lo']:.3f}", f"{f['m_hi']:.3f}",
                               str(int(f["T_edge"])), str(int(f["m_edge"])),
                               f"{sm['frac_post']:.4f}", f"{sm['frac_tract']:.4f}", str(sm["n_tracts"]),
                               f"{sm['med_tract']:.0f}", f"{sm['mean_tract']:.0f}", f"{sm['max_tract']:.0f}",
                               f"{sm['T_from_tracts']:.1f}", f"{f['ll']:.2f}",
                               f"{L.projection_q(obs[gw][:, 2 * i:2 * i + 2], M['pBe'][gw], M['pTe'][gw]):.4f}",
                               f"{sm['F_long']:.4f}", f"{np.log10(f['T']):.4f}"]))
    hdr = ["sample", "group", "site", "q26c"] + [f"T{x}_r{r}" for r in L.R_BRACKET for x in ("", "_lo", "_hi")] + \
          ["m", "m_lo", "m_hi", "T_at_grid_edge", "m_at_grid_edge", "frac_post", "frac_tract", "n_tracts",
           "median_tract_bp", "mean_tract_bp", "longest_tract_bp", "T_from_tracts_r1.4", "loglik", "q_m", "F_long",
           "logT"]
    open(f"{od}/fits.tsv", "w").write("\t".join(hdr) + "\n" + "\n".join(rows) + "\n")
    open(f"{od}/tracts.tsv", "w").write("sample\thap\tchrom\tstart\tend\tlength\n" + "\n".join(trows) + "\n")
    print(f"wrote {od}/fits.tsv", flush=True)


if __name__ == "__main__":
    main()
