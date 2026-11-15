#!/bin/bash
## 01_convert.sh -- METHODS A2: folded SweepFinder2 inputs for all five lineages x chr1..chr25.
## Site set is IDENTICAL to the unfolded run (still AA-filtered, still --polymorphic-only), so
## the folded-versus-unfolded difference is attributable to polarization alone.
set -euo pipefail
cd "$(dirname "$0")"; source common.sh
: > "${FCLR_OUT}/tasks_convert.txt"
for pop in ${POPS}; do for i in $(seq 1 25); do printf '%s chr%d\n' "$pop" "$i" >> "${FCLR_OUT}/tasks_convert.txt"; done; done
grep . "${FCLR_OUT}/tasks_convert.txt" | xargs -L 1 -P "${CONVERT_PARALLEL:-25}" bash 01_convert_worker.sh
echo "folded inputs: $(ls "${SF2DIR}"/sfs_input_folded/*.sf2.txt 2>/dev/null | wc -l) of 125"
