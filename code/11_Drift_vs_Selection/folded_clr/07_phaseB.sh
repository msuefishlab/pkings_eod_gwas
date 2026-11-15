#!/bin/bash
## 07_phaseB.sh -- Phase B: Peak 6 against the drift-made fixed-difference cluster null.
## The B1 precondition (>= 200 clusters per lineage, and the count table agreeing with the BED
## line count) is checked inside 07_phaseB.R and stops the stage if it fails.
set -euo pipefail
cd "$(dirname "$0")"; source common.sh
cd "${root}"
sx Rscript code/11_Drift_vs_Selection/folded_clr/07_phaseB.R \
   --cluster_dir "${CLUSTER_DIR}" --out_dir "${FCLR_OUT}"
