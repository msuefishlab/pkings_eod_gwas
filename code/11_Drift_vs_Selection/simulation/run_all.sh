#!/bin/bash
## run_all.sh -- Peaks 4/5 simulation null (stages 01-07): every stage, in order.
##   bash code/11_Drift_vs_Selection/simulation/run_all.sh <stage>
## Stages: setup | observed | windows | coalescent | slim-tasks | slim | all
##
## RUN = agree (default) uses the six-species-agreement polarization run of the drift test
## (fixed_differences_agree) and its seven admissible A1w models. Its params.json is also the input
## of the LD-decay scripts (23*).
## RUN = first uses the first-pass fixed_differences run and its four adequate-but-inadmissible
## A1w models (see README.md).
##
## Run the heavy stages inside one sbatch job on a large node; the fan-outs below are xargs over a
## task file rather than SLURM arrays.
set -euo pipefail
root="$(git rev-parse --show-toplevel)"
cd "${root}"
source "${root}/config.env"

dir=code/11_Drift_vs_Selection/simulation
RUN=${RUN:-agree}
case "${RUN}" in
  agree) fd_out=output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences_agree
         select=admissible; as_run=${dir}/params.as_run.agree.json
         msp_worker=${dir}/04c_worker_msprime_agree.py; gen=${dir}/05d_gen_slim_tasks_agree.py ;;
  first) fd_out=output_data/11_Drift_vs_Selection/demographic_modeling/fixed_differences
         select=adequate;   as_run=${dir}/params.as_run.json
         msp_worker=${dir}/04b_worker_msprime_positions.py; gen=${dir}/05_gen_slim_tasks.py ;;
  *) echo "ERROR: RUN must be agree or first" >&2; exit 1 ;;
esac
work=${SIM_WORK:-${sweep_scratch}/peaks45_sim/${RUN}}
out=${SIM_OUT:-output_data/11_Drift_vs_Selection/simulation}
P=${SIM_WORKERS:-128}
BUDGET=${SIM_BUDGET_S:-3600}
MSP=${MSPRIME_PYTHON:-python3}
stage=${1:-all}

setup() {
  mkdir -p "${work}"/{chunks,joint,pos,slim_out,shards,obs}
  ## The fitted histories are DERIVED from the drift-test outputs, never copied from scratch.
  python3 "${dir}/01_fitted_history.py" --agree-dir "${fd_out}" --select "${select}" \
          --out "${work}/params.json"
  ## Regression check against the as-run history, so a change in the drift test's fits is visible.
  python3 "${dir}/01_fitted_history.py" --agree-dir "${fd_out}" --select "${select}" --check "${as_run}"
}

observed() {
  ## Per-chunk scans of the same site table the drift test built. Peak-4/5 delta thresholds and the
  ## two call-rate filters; both variants are reported (METHODS.md section 4).
  ls "${fd_out}"/site_table/snps.*.tsv.gz > "${work}/obs/chunks.list"
  grep -c . "${work}/obs/chunks.list"
  xargs -a "${work}/obs/chunks.list" -P "${P}" -I{} \
    sh -c 'python3 "$0" "$1" "$2/obs/$(basename "$1" .tsv.gz).counts.tsv"' "${dir}/02_observed_counts.py" '{}' "${work}"
  xargs -a "${work}/obs/chunks.list" -P "${P}" -I{} \
    sh -c 'python3 "$0" "$1" "$2/obs/$(basename "$1" .tsv.gz).counts_callrate.tsv"' "${dir}/02b_observed_counts_callrate.py" '{}' "${work}"
}

windows() {
  ## Empirical genome-wide baseline of 354-kb / 287-kb window counts -- METHODS.md section 7 item 3.
  xargs -a "${work}/obs/chunks.list" -P "${P}" -I{} \
    sh -c 'python3 "$0" "$1" "$2/obs/$(basename "$1" .tsv.gz).windows.tsv"' "${dir}/03_window_scan.py" '{}' "${work}"
}

coalescent() {
  ## msprime null with linkage. One task per (model, mu, task id); each runs 2-Mb chunks until its
  ## wall-clock budget is spent, so coverage scales with the number of tasks, not with a chunk count.
  : > "${work}/msp_tasks.txt"
  for key in $(python3 -c 'import json,sys; print(" ".join(k for k in json.load(open(sys.argv[1])) if "|" in k))' "${work}/params.json"); do
    m=${key%%|*}; mu=${key##*|}
    for t in $(seq 1 "${SIM_MSP_TASKS:-8}"); do printf '%s\t%s\t%s\n' "${m}" "${mu}" "${t}" >> "${work}/msp_tasks.txt"; done
  done
  grep -c . "${work}/msp_tasks.txt"
  grep . "${work}/msp_tasks.txt" | xargs -P "${P}" -L 1 \
    sh -c 'exec "$0" "$1" "$2" "$3" "$4" "$5"' "${MSP}" "${root}/${msp_worker}" '{}' "${BUDGET}" "${work}"
}

slim_tasks() {
  python3 "${root}/${gen}" "${work}/tasks.txt"
  grep -c . "${work}/tasks.txt"
  split -n l/"${SIM_SHARDS:-200}" -d -a 3 "${work}/tasks.txt" "${work}/shards/shard_"
}

slim() {
  : "${SLIM_BIN:?SLIM_BIN is not set -- point it at a SLiM 5.2 binary, see README.md}"
  ls "${work}"/shards/shard_* > "${work}/shards.list"
  xargs -a "${work}/shards.list" -P "${P}" -I{} \
    sh -c 'bash "$0" "$1" "$2/slim_out/$(basename "$1").tsv" "$3" "$4"' \
      "${dir}/07b_shard_slim_parameterised.sh" '{}' "${work}" "${SLIM_BIN}" "${root}/${dir}/power.slim"
  cat "${work}"/slim_out/shard_*.tsv > "${work}/slim_power.tsv"
  wc -l "${work}/slim_power.tsv"
}

case "${stage}" in
  setup) setup ;;
  observed) observed ;;
  windows) windows ;;
  coalescent) coalescent ;;
  slim-tasks) slim_tasks ;;
  slim) slim ;;
  all) setup; observed; windows; coalescent; slim_tasks; slim ;;
  *) echo "ERROR: unknown stage ${stage}" >&2; exit 1 ;;
esac
mkdir -p "${out}"
echo "work dir: ${work}  |  results dir: ${out}"
