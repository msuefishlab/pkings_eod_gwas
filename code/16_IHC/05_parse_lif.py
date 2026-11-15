#!/usr/bin/env python3
"""Read Leica .lif files (LAS X) without readlif, and tabulate SPTBN4 IHC acquisition settings.

A .lif file is a UTF-16 XML header followed by memory blocks. Each image element names
its memory block, its X/Y/Z dimensions (with physical length) and its channels. Channel
planes sit at the `ChannelDescription BytesInc` offset inside the block. The per-channel
filter cube, exposure, lamp (IL) intensity and gain come from the `WideFieldChannelInfo`
entries of the image's own hardware-setting attachment (same order as the channels).

As a module: `read_lif(path)` returns (fh, images); `plane(fh, img, cube)` returns one
channel as a float array (cube "L5" = 488/GFP = SPTBN4, "DAP" = DAPI, "Ana" = DIC).

As a script: writes one row per image with a `comparable_group` key. Fields can be
compared on raw intensity only within one comparable_group (same session, objective,
image size, and identical L5 exposure, IL and gain).

Usage:
    source config.env
    python3 code/16_IHC/05_parse_lif.py "$sptbn4_ihc_dir" output_data/16_IHC/sptbn4_acquisition_settings.csv
"""
import csv
import os
import re
import struct
import sys
import xml.etree.ElementTree as ET

import numpy as np

# Staining/imaging sessions (2025). One APA (TP) and one BAM (BP) fish throughout;
# APA and BAM slides were stained together in each session.
SESSIONS = {
    "PKINGS_SPTBN_IHC.lif": "2025-10-14",
    "PKINGS_SPTBN_IHC_OCT23.lif": "2025-10-23",
    "PKINGS_SPTBN_IHC_NOV11_REP_SECTS_CHECKED.lif": "2025-11-11",
}
LUT_FOR_CUBE = {"L5": "Green", "DAP": "Blue", "TXR": "Red"}


def read_lif(path):
    fh = open(path, "rb")
    magic, _, _, nchar = struct.unpack("<iiBi", fh.read(13))
    if magic != 0x70:
        raise ValueError(f"{path}: not a LIF file")
    root = ET.fromstring(fh.read(nchar * 2).decode("utf-16-le"))

    blocks = {}
    while True:
        head = fh.read(9)
        if len(head) < 9:
            break
        magic, _, _ = struct.unpack("<iiB", head)
        if magic != 0x70:
            break
        size, = struct.unpack("<q", fh.read(8))
        fh.read(1)
        n, = struct.unpack("<i", fh.read(4))
        bid = fh.read(n * 2).decode("utf-16-le")
        blocks[bid] = (fh.tell(), size)
        fh.seek(size, 1)

    images = []
    for el in root.iter("Element"):
        img = el.find("Data/Image")
        mem = el.find("Memory")
        if img is None or mem is None or int(mem.get("Size", 0)) == 0:
            continue
        dims = {int(d.get("DimID")): d for d in img.find("ImageDescription/Dimensions")}
        chans = [c.attrib for c in img.find("ImageDescription/Channels")]
        # Hardware settings of this element only (not its children).
        wf = []
        for att in el.findall("Data/Image/Attachment"):
            if att.get("Name") == "HardwareSetting":
                wf = [w.attrib for w in att.iter("WideFieldChannelInfo")]
                break
        objective = next((e.get("ObjectiveName") for e in el.iter() if e.get("ObjectiveName")), "NA")
        mag = next((e.get("Magnification") for e in el.iter() if e.get("Magnification")), "NA")
        x, y = dims[1], dims[2]
        channels, lut_ok = [], True
        for i, c in enumerate(chans):
            w = wf[i] if i < len(wf) else {}
            cube = w.get("FluoCubeName", "NA")
            if cube in LUT_FOR_CUBE and c.get("LUTName") != LUT_FOR_CUBE[cube]:
                lut_ok = False  # channel/hardware-setting order disagree (seen only in processed images)
            channels.append(dict(
                cube=cube, method=w.get("ContrastingMethodName", "NA"),
                exposure_s=float(w.get("ExposureTime", "nan")),
                IL=w.get("IL_Light-Intensity", "NA"), gain=w.get("Gain", "NA"),
                bytes_inc=int(c["BytesInc"]), bits=int(c["Resolution"])))
        images.append(dict(
            name=el.get("Name"), objective=" ".join(objective.split()), magnification=mag,
            nx=int(x.get("NumberOfElements")), ny=int(y.get("NumberOfElements")),
            nz=int(dims[3].get("NumberOfElements")) if 3 in dims else 1,
            um_per_px=float(x.get("Length")) * 1e6 / int(x.get("NumberOfElements")),
            block=blocks[mem.get("MemoryBlockID")], channels=channels, lut_ok=lut_ok))
    return fh, images


def plane(fh, img, cube):
    """One channel (first z-plane) as float64, shape (ny, nx)."""
    if not img["lut_ok"]:
        raise ValueError(f"{img['name']}: channel order does not match the hardware setting")
    ch = next(c for c in img["channels"] if c["cube"] == cube)
    if ch["bits"] != 8:
        raise ValueError(f"{img['name']}: {ch['bits']}-bit channel not supported")
    offset, _ = img["block"]
    fh.seek(offset + ch["bytes_inc"])
    a = np.frombuffer(fh.read(img["nx"] * img["ny"]), dtype=np.uint8)
    return a.reshape(img["ny"], img["nx"]).astype(np.float64)


def parse_name(name):
    """Population and control status from JRG's field names, e.g. APA_EO_SLIDE8_SECT3_SPTBN4_POS_40x,
    BAMBOMO_SPTBN_NEG_20x, APA_EO_NO_PRIMARY, 100X_BAM_POS_1. Unnamed/processed images get NA."""
    n = name.upper()
    tokens = re.split(r"[_\- ]", n)
    pop = "APA" if "APA" in tokens else "BAM" if ("BAM" in tokens or "BAMBOMO" in tokens) else "NA"
    if any(re.fullmatch(r"NEG\d*", t) for t in tokens) or "NO_PRIMARY" in n:
        control = "neg"
    elif any(re.fullmatch(r"POS\d*", t) for t in tokens) or "PRIMARY" in tokens:
        control = "pos"
    else:
        control = "NA"
    return pop, control


def settings_rows(lif_dir):
    rows = []
    for fname, session in SESSIONS.items():
        fh, images = read_lif(os.path.join(lif_dir, fname))
        for img in images:
            pop, control = parse_name(img["name"])
            l5 = next((c for c in img["channels"] if c["cube"] == "L5"), None)
            dapi = next((c for c in img["channels"] if c["cube"] == "DAP"), None)
            processed = img["name"].startswith(("Multifocus", "Depth Map"))
            usable = (pop != "NA" and control != "NA" and img["nz"] == 1 and not processed
                      and img["lut_ok"] and l5 is not None and dapi is not None)
            group = "NA"
            if usable:
                group = (f"{session}|{img['magnification']}x|{img['nx']}px|"
                         f"L5 {l5['exposure_s']:g}s IL{l5['IL']} g{l5['gain']}")
            rows.append(dict(
                lif=fname, session=session, field=img["name"], population=pop, control=control,
                objective=img["objective"], magnification=img["magnification"], nx=img["nx"],
                nz=img["nz"], lut_ok=img["lut_ok"], um_per_px=f"{img['um_per_px']:.5f}",
                l5_exposure_s=l5["exposure_s"] if l5 else "NA", l5_IL=l5["IL"] if l5 else "NA",
                l5_gain=l5["gain"] if l5 else "NA",
                dapi_exposure_s=dapi["exposure_s"] if dapi else "NA",
                usable=usable, comparable_group=group))
        fh.close()
    return rows


def main(lif_dir, out_csv):
    rows = settings_rows(lif_dir)
    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} images ({sum(r['usable'] for r in rows)} usable fields) to {out_csv}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
