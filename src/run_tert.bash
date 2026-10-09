#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# run_tert.bash
#
# Search each species in inputs/species.tsv (src/select_species.py) for
# TERT, as one LSF array (src/worker_tert.bash). Species with an existing
# outputs/tert/<species>.rbd.domtbl are skipped (unless FORCE=1).
# ----------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

SPECIES="${SPECIES:-${ROOT}/inputs/species.tsv}"
export MINIPROT="${MINIPROT:-${ROOT}/bin/miniprot}"
export HMMSEARCH="${HMMSEARCH:-/software/team301/hmmer-3.4/src/hmmsearch}"
export PANEL="${PANEL:-${ROOT}/inputs/tert_panel.faa}"
export RBD_HMM="${RBD_HMM:-${ROOT}/inputs/Telomerase_RBD.hmm}"
export OUT_DIR="${OUT_DIR:-${ROOT}/outputs/tert}"
export THREADS="${THREADS:-8}"
LOGDIR="${ROOT}/outputs/logs/tert"
mkdir -p "$OUT_DIR" "$LOGDIR"

memGB="${memGB:-16}"
QUEUE="${QUEUE:-normal}"
memMB=$(( memGB * 1024 ))

[[ -x "$MINIPROT" && -x "$HMMSEARCH" ]] || { echo "[run_tert] missing miniprot or hmmsearch" >&2; exit 1; }

# the species still to do, with a header so worker index i is line i + 1
todo="${LOGDIR}/todo.$(date +%Y%m%d%H%M%S).tsv"
{
    head -1 "$SPECIES"
    tail -n +2 "$SPECIES" | while IFS=$'\t' read -r species fasta rest; do
        if [[ -s "${OUT_DIR}/${species}.rbd.domtbl" && "${FORCE:-0}" != "1" ]]; then
            continue
        fi
        printf '%s\t%s\t%s\n' "$species" "$fasta" "$rest"
    done
} > "$todo"
n=$(( $(wc -l < "$todo") - 1 ))
if (( n == 0 )); then
    echo "[run_tert] nothing to do"
    rm "$todo"
    exit 0
fi

bsub \
  -J "tert[1-${n}]%200" \
  -q "$QUEUE" \
  -n "$THREADS" \
  -R "span[hosts=1] select[mem>${memMB}] rusage[mem=${memMB}]" \
  -M "$memMB" \
  -o "${LOGDIR}/%J.%I.out" \
  -e "${LOGDIR}/%J.%I.err" \
  bash "${SCRIPT_DIR}/worker_tert.bash" "$todo"

echo "[run_tert] Submitted an array of ${n} species (todo list: ${todo})."
