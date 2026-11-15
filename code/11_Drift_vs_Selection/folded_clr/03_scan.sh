#!/bin/bash
## 03_scan.sh -- genome-wide folded CLR scan, five lineages x chr1..chr25 (METHODS A4).
set -euo pipefail
cd "$(dirname "$0")"; source common.sh
: > "${FCLR_OUT}/tasks_scan.txt"
for pop in ${POPS}; do for i in $(seq 1 25); do printf '%s chr%d\n' "$pop" "$i" >> "${FCLR_OUT}/tasks_scan.txt"; done; done
grep . "${FCLR_OUT}/tasks_scan.txt" | xargs -L 1 -P "${SCAN_PARALLEL:-25}" bash 03_scan_worker.sh
echo "folded CLR scans: $(ls "${SF2DIR}"/clr_folded/*.clr.txt 2>/dev/null | wc -l) of 125"
printf 'lineage\tn_gridpoints\n' > "${FCLR_OUT}/folded_scan_gridpoints.tsv"
for pop in ${POPS}; do
  n=$(cat "${SF2DIR}"/clr_folded/${pop}.chr*.clr.txt 2>/dev/null | grep -vc '^location' || true)
  printf '%s\t%s\n' "$pop" "$n" >> "${FCLR_OUT}/folded_scan_gridpoints.tsv"
done
cat "${FCLR_OUT}/folded_scan_gridpoints.tsv"
