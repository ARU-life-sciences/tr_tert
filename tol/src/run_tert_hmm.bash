#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# run_tert_hmm.bash
#
# The sensitive TERT search (profile HMM against six-frame ORFs), as one
# LSF array (src/worker_tert_hmm.bash), for the species in SPECIES
# (species<TAB>fasta..., with a header). Species with an existing
# outputs/tert_hmm/<species>.tert.domtbl are skipped (unless FORCE=1).
# ----------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

SPECIES="${SPECIES:-${ROOT}/inputs/species.tsv}"
export ESL_TRANSLATE="${ESL_TRANSLATE:-/software/team301/hmmer-3.4/easel/miniapps/esl-translate}"
export TERT_HMM="${TERT_HMM:-${ROOT}/inputs/tert_round3.hmm}"
export HMMSEARCH="${HMMSEARCH:-/software/team301/hmmer-3.4/src/hmmsearch}"
export RBD_HMM="${RBD_HMM:-${ROOT}/inputs/Telomerase_RBD.hmm}"
export OUT_DIR="${OUT_DIR:-${ROOT}/outputs/tert_hmm}"
export THREADS="${THREADS:-8}"
LOGDIR="${ROOT}/outputs/logs/tert_hmm"
mkdir -p "$OUT_DIR" "$LOGDIR"

memGB="${memGB:-16}"
QUEUE="${QUEUE:-normal}"
memMB=$(( memGB * 1024 ))

[[ -x "$ESL_TRANSLATE" && -x "$HMMSEARCH" && -s "$TERT_HMM" ]] || { echo "[run_tert_hmm] missing esl-translate, hmmsearch or the TERT HMM" >&2; exit 1; }

# the species still to do, with a header so worker index i is line i + 1
todo="${LOGDIR}/todo.$(date +%Y%m%d%H%M%S).tsv"
{
    head -1 "$SPECIES"
    tail -n +2 "$SPECIES" | while IFS=$'\t' read -r species fasta rest; do
        if [[ -s "${OUT_DIR}/${species}.tert.domtbl" && "${FORCE:-0}" != "1" ]]; then
            continue
        fi
        printf '%s\t%s\t%s\n' "$species" "$fasta" "$rest"
    done
} > "$todo"
n=$(( $(wc -l < "$todo") - 1 ))
if (( n == 0 )); then
    echo "[run_tert_hmm] nothing to do"
    rm "$todo"
    exit 0
fi

bsub \
  -J "terthmm[1-${n}]%200" \
  -q "$QUEUE" \
  -n "$THREADS" \
  -R "span[hosts=1] select[mem>${memMB}] rusage[mem=${memMB}]" \
  -M "$memMB" \
  -o "${LOGDIR}/%J.%I.out" \
  -e "${LOGDIR}/%J.%I.err" \
  bash "${SCRIPT_DIR}/worker_tert_hmm.bash" "$todo"

echo "[run_tert_hmm] Submitted an array of ${n} species (todo list: ${todo})."
