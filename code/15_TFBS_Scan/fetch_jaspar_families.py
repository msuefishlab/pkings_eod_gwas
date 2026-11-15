#!/usr/bin/env python3
"""
Fetch the JASPAR TF class and family for every motif in the MEME motif file.

Writes input_data/15_TFBS_Scan/jaspar2024_families.tsv with columns
motif_id, motif_name, tf_class, tf_family. Multiple families are joined with ';'.
Empty values become 'Unknown'. Run once; the table is the cached copy used by
matched_null.py and candidate_variants.py.

Usage:
  python3 code/15_TFBS_Scan/fetch_jaspar_families.py \\
    --motifs input_data/15_TFBS_Scan/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt \\
    --out input_data/15_TFBS_Scan/jaspar2024_families.tsv
"""

import argparse
import csv
import json
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

API = "https://jaspar.elixir.no/api/v1/matrix/{}/?format=json"


def parse_args():
    p = argparse.ArgumentParser(description="Fetch JASPAR TF class/family per motif.")
    p.add_argument("--motifs", required=True, help="MEME-format JASPAR motif file")
    p.add_argument("--out", required=True, help="Output TSV")
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


def read_motifs(path):
    motifs = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("MOTIF"):
                parts = line.split()
                motifs.append((parts[1], parts[2] if len(parts) > 2 else parts[1]))
    return motifs


def fetch(motif_id, retries=4):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(API.format(motif_id), timeout=30) as r:
                d = json.load(r)
            # JASPAR strings sometimes carry stray '\r' or spaces; strip them
            fam = ";".join(s.strip() for s in d.get("family") or [] if s.strip()) or "Unknown"
            cls = ";".join(s.strip() for s in d.get("class") or [] if s.strip()) or "Unknown"
            return cls, fam
        except Exception as e:  # network hiccup: back off and retry
            if attempt == retries - 1:
                print(f"[ERROR] {motif_id}: {e}", file=sys.stderr)
                return "Unknown", "Unknown"
            time.sleep(2 ** attempt)


def main():
    args = parse_args()
    motifs = read_motifs(args.motifs)
    with ThreadPoolExecutor(args.threads) as ex:
        results = list(ex.map(lambda m: fetch(m[0]), motifs))
    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["motif_id", "motif_name", "tf_class", "tf_family"])
        for (mid, name), (cls, fam) in zip(motifs, results):
            w.writerow([mid, name, cls, fam])
    n_unknown = sum(1 for _, fam in results if fam == "Unknown")
    print(f"[INFO] Wrote {len(motifs)} motifs to {args.out} ({n_unknown} with Unknown family)")


if __name__ == "__main__":
    main()
