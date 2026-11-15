#!/bin/bash
## run_fimo_units.sh
## FIMO scan of all unit haplotype sequences (peak + matched-null pool), in parallel chunks.
## See code/15_TFBS_Scan/README.md, stage F. Runs locally; needs MEME Suite 5.5.8 (fimo, fasta-get-markov)
## and a python3 with numpy and pandas. Set MEME_BIN and PY if they are not on the PATH.
##
## --no-pgc stops FIMO parsing coordinates out of the headers, which would offset the hits;
## --text streams hits without the 100k stored-score cap and without q-values.
##
## Usage:
##   bash code/15_TFBS_Scan/run_fimo_units.sh [N_CHUNKS]

set -eo pipefail

root="$(git rev-parse --show-toplevel)"
source "${root}/config.env"
set -u

N_CHUNKS="${1:-14}"
OUT="${root}/output_data/15_TFBS_Scan"
CHUNK_DIR="${OUT}/fimo_chunks"
MEME_BIN="${MEME_BIN:-$(dirname "$(command -v fimo)")}"
PY="${PY:-$(command -v python3)}"
PARSER="${root}/code/15_TFBS_Scan/parse_fimo_hits.py"

[[ -s "${OUT}/genome.bg" ]] || "${MEME_BIN}/fasta-get-markov" -m 0 -dna "${TFBS_REF}" "${OUT}/genome.bg"

rm -rf "${CHUNK_DIR}"
mkdir -p "${CHUNK_DIR}"

# Split units.fa into N_CHUNKS files by record count (records are 2 lines: header + sequence)
awk -v n="${N_CHUNKS}" -v dir="${CHUNK_DIR}" '
  /^>/ { i++; f = sprintf("%s/chunk_%02d.fa", dir, i % n) }
  { print > f }' "${OUT}/units.fa"

export JASPAR_MOTIFS OUT MEME_BIN PY PARSER
# BSD xargs caps -I replacements at 255 bytes, so pass the chunk path as $1 instead
printf '%s\n' "${CHUNK_DIR}"/chunk_*.fa | xargs -P "${N_CHUNKS}" -n 1 bash -c '
  set -eo pipefail
  fa="$1"; base="${fa%.fa}"
  "${MEME_BIN}/fimo" --text --no-pgc --thresh 1e-3 --bfile "${OUT}/genome.bg" \
      "${JASPAR_MOTIFS}" "${fa}" 2> "${base}.fimo.log" \
    | "${PY}" "${PARSER}" parse --units "${OUT}/units.tsv" --fasta "${fa}" \
        --motifs "${JASPAR_MOTIFS}" --out "${base}.npz" 2> "${base}.parse.log"
' _

"${PY}" "${PARSER}" merge --motifs "${JASPAR_MOTIFS}" --out "${OUT}/hits.npz" "${CHUNK_DIR}"/chunk_*.npz
cat "${CHUNK_DIR}"/*.parse.log
