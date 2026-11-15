#!/usr/bin/env python3
"""Image masks used by 06_quantify_sptbn4_ihc.py.

band_mask: the DAPI-dense band (rows of electrocyte nuclei, grown by BAND_DILATE_UM to the electrocyte faces),
built from DAPI only, so it does not depend on the antibody signal being measured.
blob_mask: compact bright objects (vessel cross-sections, cells) in the background-subtracted antibody channel.
"""
import numpy as np
from scipy import ndimage as ndi

BAND_SIGMA_UM = 7.0                # smoothing of the DAPI density map (5-10 um give the same answer)
BAND_DILATE_UM = 4.0               # grow the nuclear band to include the electrocyte faces
BLOB_MAX_AREA_UM2 = 300.0          # compact bright objects (vessel/RBC-sized)
BLOB_MAX_EXTENT_UM = 30.0


def otsu(values, bins=512):
    hist, edges = np.histogram(values, bins=bins)
    mids = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * mids) / np.maximum(w0, 1)
    m1 = (np.sum(hist * mids) - np.cumsum(hist * mids)) / np.maximum(w1, 1)
    return mids[np.argmax(w0 * w1 * (m0 - m1) ** 2)]


def band_mask(dapi, um_per_px):
    """DAPI-dense band: smoothed nuclear density above its Otsu threshold."""
    # Clip bright DAPI patches (debris, autofluorescence) so they do not set the threshold.
    d = np.clip(dapi, *np.percentile(dapi, [5, 99]))
    dens = ndi.gaussian_filter(d, BAND_SIGMA_UM / um_per_px)
    band = dens > otsu(dens.ravel())
    band = ndi.binary_opening(band, iterations=2)
    return ndi.binary_dilation(band, iterations=max(1, int(round(BAND_DILATE_UM / um_per_px))))


def blob_mask(bs, thr, um_per_px):
    """Compact bright components. Thresholding at the bright-structure level keeps the
    membrane bands long and connected, so only isolated objects (vessel cross-sections,
    cells) pass the size/extent filter."""
    lab, n = ndi.label(bs > thr)
    if n == 0:
        return np.zeros_like(bs, dtype=bool)
    keep = np.zeros(n + 1, dtype=bool)
    max_area = BLOB_MAX_AREA_UM2 / um_per_px ** 2
    max_ext = BLOB_MAX_EXTENT_UM / um_per_px
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        area = np.sum(lab[sl] == i)
        ext = max(sl[0].stop - sl[0].start, sl[1].stop - sl[1].start)
        keep[i] = area <= max_area and ext <= max_ext
    return ndi.binary_dilation(keep[lab], iterations=2)
