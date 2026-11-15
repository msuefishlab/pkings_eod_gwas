#!/bin/bash
## common.sh -- shared environment for the Peak 6 folded-CLR and fixed-difference-null work.
## Sourced by every numbered stage. Every path below derives from the repository root; the
## singularity bind (SINGULARITY_BIND) defaults to the scratch directory and is overridable.
set -uo pipefail
root=${root:-$(git rev-parse --show-toplevel)}
source "${root}/config.env"
export root

export SINGULARITY_BIND="${SINGULARITY_BIND:-${scratch_dir}}"
sx()  { singularity exec "${gwas_tools_image}" "$@"; }
sf2() { singularity exec "${sweepfinder2_image}" "$@"; }

## read-only inputs produced by the sweep-detection pipeline
export SD=${SD:-${root}/output_data/11_Drift_vs_Selection/sweep_detection}
export EV=${EV:-${SD}/evidence}
export SF2DIR=${SF2DIR:-${SD}/sweepfinder2}
export CLUSTER_DIR=${CLUSTER_DIR:-${root}/output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree}

## this analysis's own outputs
export FCLR_OUT=${FCLR_OUT:-${root}/output_data/11_Drift_vs_Selection/folded_clr}
export LOGDIR=${LOGDIR:-${FCLR_OUT}/logs}
export POPS=${POPS:-"BP1 TP1 BP2 TP2 BP3"}
mkdir -p "${FCLR_OUT}" "${LOGDIR}" "${FCLR_OUT}/audit" \
         "${SF2DIR}/sfs_input_folded" "${SF2DIR}/background_folded" "${SF2DIR}/clr_folded"
