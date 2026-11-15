#!/bin/bash
## 04_merge_folded.sh -- merge the per-chromosome folded scans into
## sweepfinder2/merged/<pop>.genome_wide_clr_folded.txt, which is the filename
## region_stats_lib.R's clr_folded registry entry already expects.
##
## This is the repository's own 13_sf2_merge.sh in its `folded` mode -- the column layout is
## therefore identical to the unfolded merge by construction, which is what METHODS A4 requires.
set -euo pipefail
cd "$(dirname "$0")"; source common.sh
bash "${root}/code/11_Drift_vs_Selection/sweep_detection/sweepfinder2/13_sf2_merge.sh" folded
