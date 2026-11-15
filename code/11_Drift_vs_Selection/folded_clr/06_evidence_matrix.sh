#!/bin/bash
## 06_evidence_matrix.sh -- the decision rule applied twice over the same region statistics:
##   (a) CLR_STAT=clr        -> the re-derived unfolded run, which is the self-check of A5;
##   (b) CLR_STAT=clr_folded -> the folded primary triple of A6.
## Both write into output_data/11_Drift_vs_Selection/folded_clr, never over the headline evidence tables.
set -euo pipefail
cd "$(dirname "$0")"; source common.sh
cd "${root}"
S=code/11_Drift_vs_Selection/folded_clr/06_evidence_matrix_clrstat.R

CLR_STAT=clr        sx Rscript "$S" --suffix _folded --out_suffix _rederived_unfolded --out_dir "${FCLR_OUT}"
CLR_STAT=clr_folded sx Rscript "$S" --suffix _folded --out_suffix _folded             --out_dir "${FCLR_OUT}"
ls -lh "${FCLR_OUT}"/evidence_matrix_*.tsv "${FCLR_OUT}"/calibration_*.tsv
