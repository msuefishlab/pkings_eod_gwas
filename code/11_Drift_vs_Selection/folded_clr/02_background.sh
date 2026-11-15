#!/bin/bash
## 02_background.sh -- per-lineage empirical folded background spectrum (METHODS A4).
set -euo pipefail
cd "$(dirname "$0")"; source common.sh
for pop in ${POPS}; do bash 02_background_worker.sh "$pop"; done
