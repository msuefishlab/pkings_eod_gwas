#!/usr/bin/env python3
"""Per-field SPTBN4 immunofluorescence metrics from the 2025 Leica .lif files (APA TP vs BAM BP; Fig. 4B).

Measures the L5 (488 nm; anti-SPTBN4 + Alexa 488) signal in raw 8-bit camera images, comparing APA
and BAM only within a comparable_group (same session, objective, image size and identical L5 exposure,
lamp and gain; 05_parse_lif.py).

Masks (ihc_masks.py):
  band   DAPI-dense band (rows of electrocyte nuclei, grown by BAND_DILATE_UM to the
         electrocyte faces). Built from DAPI only, so independent of the antibody signal.
  blobs  compact bright L5 objects. In APA these are mostly stalks, plus vessels and cells, so
         blob metrics are descriptive only (not a vessel measure); defined from the antibody signal.

Primary metric: band signal = mean L5 inside the band minus the mean L5 inside the band of the
same population's no-primary control in the same group. Effect size per group: BAM/APA ratio of field medians, with a bootstrap interval
over fields. One fish per population, so fields are technical replicates: no population test.

Background sensitivity: band_signal_commonbg subtracts the mean of the APA and BAM controls'
band means from both populations; band_raw subtracts nothing. On 2025-10-23 the APA control
is brighter than the BAM control (median 45 vs 17), so that group depends on this choice.

Groups without a no-primary control (2025-10-14, 100x) are reported on raw intensity only.

Usage:
    source config.env
    python3 code/16_IHC/06_quantify_sptbn4_ihc.py "$sptbn4_ihc_dir" \
        output_data/16_IHC/sptbn4_acquisition_settings.csv output_data/16_IHC
"""
import collections
import csv
import importlib.util
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
N_BOOT = 10000
SATURATED = 255


def load(name, fname):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


ihc = load("ihc_masks", "ihc_masks.py")
lif = load("parse_lif", "05_parse_lif.py")


def band_at(dapi, um_per_px, dilate_um):
    old = ihc.BAND_DILATE_UM
    ihc.BAND_DILATE_UM = dilate_um
    try:
        return ihc.band_mask(dapi, um_per_px)
    finally:
        ihc.BAND_DILATE_UM = old


def save_qc(path, dapi, l5_bs, band, blobs, disp_max):
    """DAPI (auto-scaled) | L5 above background at the group's fixed scale | band outline
    (cyan) and vessel/cell blobs (magenta) on the L5 panel."""
    d = np.clip(dapi / max(np.percentile(dapi, 99.5), 1), 0, 1)
    g = (np.clip(l5_bs / disp_max, 0, 1) * 255).astype(np.uint8)
    edge = band ^ ndi.binary_erosion(band, iterations=1)
    dp = np.stack([np.zeros_like(g), np.zeros_like(g), (d * 255).astype(np.uint8)], axis=-1)
    dp[edge] = (0, 220, 255)
    lp = np.stack([g, g, g], axis=-1)
    ov = lp.copy()
    ov[edge] = (0, 220, 255)
    ov[blobs & ~edge] = (255, 0, 200)
    sep = np.full((g.shape[0], 6, 3), 255, np.uint8)
    Image.fromarray(np.hstack([dp, sep, lp, sep, ov])).save(path)


def boot_ratio(apa, bam, rng):
    apa, bam = np.asarray(apa), np.asarray(bam)
    ra = np.median(rng.choice(apa, (N_BOOT, apa.size)), axis=1)
    rb = np.median(rng.choice(bam, (N_BOOT, bam.size)), axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        r = rb / ra
    return np.percentile(r[np.isfinite(r)], [2.5, 97.5])


def main(lif_dir, settings_csv, out_dir):
    rng = np.random.default_rng(20261114)
    os.makedirs(os.path.join(out_dir, "qc"), exist_ok=True)
    settings = [r for r in csv.DictReader(open(settings_csv)) if r["usable"] == "True"]
    groups = collections.defaultdict(list)
    for r in settings:
        groups[r["comparable_group"]].append(r)

    handles = {}

    def get(r):
        if r["lif"] not in handles:
            fh, images = lif.read_lif(os.path.join(lif_dir, r["lif"]))
            handles[r["lif"]] = (fh, {im["name"]: im for im in images})
        fh, images = handles[r["lif"]]
        im = images[r["field"]]
        return lif.plane(fh, im, "L5"), lif.plane(fh, im, "DAP"), im["um_per_px"]

    field_rows, group_rows = [], []
    for gkey, rows in groups.items():
        pos = {p: [r for r in rows if r["population"] == p and r["control"] == "pos"] for p in ("APA", "BAM")}
        neg = {p: [r for r in rows if r["population"] == p and r["control"] == "neg"] for p in ("APA", "BAM")}
        if not (pos["APA"] and pos["BAM"]):
            continue  # one population only: nothing to compare
        has_ctrl = all(len(neg[p]) == 1 for p in ("APA", "BAM"))
        # Settings check: every field in the group shares the L5 exposure, lamp and gain.
        assert len({(r["l5_exposure_s"], r["l5_IL"], r["l5_gain"], r["nx"]) for r in rows}) == 1, gkey

        # Background per population from its no-primary control (whole field and inside the band).
        bg = {}
        for p in ("APA", "BAM"):
            if has_ctrl:
                l5n, dn, um = get(neg[p][0])
                bg[p] = dict(field=float(np.median(l5n)), thr=float(np.percentile(l5n, 99.5)),
                             band={dil: float(l5n[band_at(dn, um, dil)].mean()) for dil in (2.0, 4.0, 6.0)},
                             neg_field=neg[p][0]["field"])
            else:
                bg[p] = dict(field=0.0, thr=np.nan, band={2.0: 0.0, 4.0: 0.0, 6.0: 0.0}, neg_field="none")

        stained = [(r, *get(r)) for p in ("APA", "BAM") for r in pos[p]]
        pooled = np.concatenate([(l5 - bg[r["population"]]["field"]).ravel() for r, l5, _, _ in stained])
        bright_thr = float(np.percentile(pooled, 95))
        disp_max = float(np.percentile(pooled, 99.5))

        for r, l5, dapi, um in stained:
            b = bg[r["population"]]
            bs = l5 - b["field"]
            band = band_at(dapi, um, 4.0)
            blobs = ihc.blob_mask(bs, bright_thr, um)
            sat = l5 >= SATURATED
            common_bg = float(np.mean([bg[p]["band"][4.0] for p in ("APA", "BAM")]))
            row = dict(group=gkey, session=r["session"], field=r["field"], population=r["population"],
                       has_control=has_ctrl, neg_control=b["neg_field"], um_per_px=um,
                       saturated_frac=float(sat.mean()), band_area_frac=float(band.mean()),
                       blob_area_frac=float(blobs.mean()), bg_field=b["field"], bg_band=b["band"][4.0],
                       band_signal=float(l5[band].mean() - b["band"][4.0]),
                       band_raw=float(l5[band].mean()),
                       band_signal_commonbg=float(l5[band].mean() - common_bg),
                       band_signal_dil2=float(l5[band_at(dapi, um, 2.0)].mean() - b["band"][2.0]),
                       band_signal_dil6=float(l5[band_at(dapi, um, 6.0)].mean() - b["band"][6.0]),
                       band_signal_nosat=float(l5[band & ~sat].mean() - b["band"][4.0]),
                       band_signal_noblob=float(l5[band & ~blobs].mean() - b["band"][4.0]),
                       field_signal=float(l5.mean() - b["field"]),
                       interband_signal=float(bs[~band].mean()),
                       pos_frac=float((l5 > b["thr"]).mean()) if has_ctrl else np.nan,
                       band_pos_frac=float((l5[band] > b["thr"]).mean()) if has_ctrl else np.nan,
                       blob_mean=float(bs[blobs].mean()) if blobs.any() else np.nan,
                       p99=float(np.percentile(l5, 99)))
            field_rows.append(row)
            save_qc(os.path.join(out_dir, "qc", f"sptbn4_{r['session']}_{r['field']}.png"),
                    dapi, bs, band, blobs, disp_max)

        mets = ["band_signal", "band_raw", "band_signal_commonbg", "band_signal_dil2", "band_signal_dil6", "band_signal_nosat",
                "band_signal_noblob", "field_signal", "interband_signal", "pos_frac", "band_pos_frac",
                "blob_mean", "p99"]
        g = dict(group=gkey, has_control=has_ctrl, n_apa=len(pos["APA"]), n_bam=len(pos["BAM"]),
                 bright_thr=bright_thr, display_max=disp_max)
        for m in mets:
            va = [f[m] for f in field_rows if f["group"] == gkey and f["population"] == "APA"]
            vb = [f[m] for f in field_rows if f["group"] == gkey and f["population"] == "BAM"]
            if np.all(np.isnan(va)):
                g[f"{m}_APA"] = g[f"{m}_BAM"] = g[f"{m}_ratio"] = np.nan
                continue
            ma, mb = float(np.nanmedian(va)), float(np.nanmedian(vb))
            g[f"{m}_APA"], g[f"{m}_BAM"] = ma, mb
            g[f"{m}_ratio"] = mb / ma if ma else np.nan
            if m == "band_signal":
                lo, hi = boot_ratio(va, vb, rng)
                g["band_signal_ratio_lo"], g["band_signal_ratio_hi"] = float(lo), float(hi)
                g["band_signal_BAM_max_over_APA_min"] = max(vb) / min(va) if min(va) else np.nan
        group_rows.append(g)

    for name, rows in (("sptbn4_field_metrics.csv", field_rows), ("sptbn4_group_summary.csv", group_rows)):
        with open(os.path.join(out_dir, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            for row in rows:
                w.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in row.items()})

    print(f"{'group':40s} nA nB  band APA  band BAM  ratio [95% boot]   raw  common  dil2  dil6  nosat  field  pos_frac")
    for g in group_rows:
        print(f"{g['group']:40s} {g['n_apa']:2d} {g['n_bam']:2d} {g['band_signal_APA']:9.2f} {g['band_signal_BAM']:9.2f}"
              f"  {g['band_signal_ratio']:.2f} [{g['band_signal_ratio_lo']:.2f}, {g['band_signal_ratio_hi']:.2f}]"
              f"  {g['band_raw_ratio']:.2f}  {g['band_signal_commonbg_ratio']:.2f}  {g['band_signal_dil2_ratio']:.2f}"
              f"  {g['band_signal_dil6_ratio']:.2f}  {g['band_signal_nosat_ratio']:.2f}  {g['field_signal_ratio']:.2f}"
              f"  {g['pos_frac_ratio']:.2f}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
