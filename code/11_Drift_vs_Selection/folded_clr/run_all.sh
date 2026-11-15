#!/bin/bash
## run_all.sh -- Peak 6 folded-spectrum CLR (part A) and the fixed-difference null (part B).
## Usage (from anywhere in the repository):
##   bash code/11_Drift_vs_Selection/folded_clr/run_all.sh
## Stages 1 and 3 are the long poles (125 conversions and 125 SweepFinder2 scans); run on a
## compute node. Stage 0 is blocking: if the folded conversion is not polarization-free, nothing
## downstream runs.
set -euo pipefail
cd "$(dirname "$0")"; source common.sh
echo "root=${root}"; echo "out =${FCLR_OUT}"

bash 00_folded_precondition.sh BP2 chr1
bash 01_convert.sh
bash 02_background.sh
bash 03_scan.sh
bash 04_merge_folded.sh
bash 05_region_stats.sh
bash 06_evidence_matrix.sh
bash 07_phaseB.sh
bash 08_correction.sh
ls -lh "${FCLR_OUT}"
