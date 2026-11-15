#!/bin/bash
## 08_correction.sh -- folded beside unfolded for all six peaks, the blocking reproduction
## check, and the six-test (Bonferroni and Benjamini-Hochberg) correction.
set -euo pipefail
cd "$(dirname "$0")"; source common.sh
cd "${root}"
sx Rscript code/11_Drift_vs_Selection/folded_clr/08_correction.R \
   --out_dir "${FCLR_OUT}" \
   --rederived "${FCLR_OUT}/evidence_matrix_rederived_unfolded.tsv" \
   --folded    "${FCLR_OUT}/evidence_matrix_folded.tsv" \
   --published "${EV}/evidence_matrix.tsv"
