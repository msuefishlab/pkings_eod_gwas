#!/bin/bash
## 05_region_stats.sh -- METHODS A5: recalibrate with the IDENTICAL headline operator.
## 02_region_stats.R is the repository's own script, run unchanged, with clr_folded added to the
## statistic list. It writes region_stats_folded.tsv.gz and region_stats_operators_folded.tsv
## into the sweep-detection evidence directory, beside the other sensitivity rows (_n16,
## _r50, _rsbgw, _w40), which is where that script's --suffix outputs belong.
set -euo pipefail
cd "$(dirname "$0")"; source common.sh
cd "${root}"
sx Rscript code/11_Drift_vs_Selection/sweep_detection/evidence/02_region_stats.R \
   --stats rsb,clr,pi_ratio,clr_folded --suffix _folded
cp -f "${EV}/region_stats_operators_folded.tsv" "${FCLR_OUT}/region_stats_operators.tsv"
