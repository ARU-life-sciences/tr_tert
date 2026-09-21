#!/usr/bin/env bash
set -euo pipefail

# ./run_confirm_positional_clade.bash <clade_name> [top_n] [window_bp]
#
# Submits one LSF job per species running confirm_repeat_positional_clade.py
# (top-N Stage-0 candidates each, via `tidk search`, restricted to
# chromosome termini) - see that script's docstring for why this is a
# separate, necessary step from Stage 0 (`tidk explore`'s top hit is not
# necessarily telomeric). Skips species that already have a summary file
# (unless FORCE=1).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLADE_ROOT="$(dirname "$SCRIPT_DIR")"

CLADE_NAME="${1:?usage: run_confirm_positional_clade.bash <clade_name> [top_n] [window_bp]}"
TOP_N="${2:-3}"
WINDOW="${3:-5000}"

ASSEMBLY_LIST="${CLADE_ROOT}/inputs/${CLADE_NAME}_paths.txt"
[[ -s "$ASSEMBLY_LIST" ]] || { echo "No assembly list at $ASSEMBLY_LIST" >&2; exit 1; }

SUMMARIES_DIR="${CLADE_ROOT}/outputs/positional/summaries"
LOGDIR="${CLADE_ROOT}/outputs/logs/positional"
mkdir -p "$SUMMARIES_DIR" "$LOGDIR"

QUEUE="${QUEUE:-normal}"
memGB="${memGB:-4}"
nCPU="${nCPU:-1}"
memMB=$(( memGB * 1024 ))
FORCE="${FORCE:-0}"

echo "[run_confirm_positional_clade] clade: $CLADE_NAME, top_n=$TOP_N, window=$WINDOW"

while IFS=$'\t' read -r species fasta; do
  [[ -z "$species" || -z "$fasta" ]] && continue
  summary="${SUMMARIES_DIR}/${species}.tsv"
  if [[ "$FORCE" != "1" && -s "$summary" ]]; then
    echo "[run_confirm_positional_clade] ${species}: exists, skipping."
    continue
  fi

  log_out="${LOGDIR}/${species}.%J.out"
  log_err="${LOGDIR}/${species}.%J.err"

  bsub \
    -J "confirmpos_${CLADE_NAME}_${species}" \
    -q "$QUEUE" \
    -n "$nCPU" \
    -R "span[hosts=1] select[mem>${memMB}] rusage[mem=${memMB}]" \
    -M "$memMB" \
    -o "$log_out" \
    -e "$log_err" \
    python3 "${SCRIPT_DIR}/confirm_repeat_positional_clade.py" "$CLADE_NAME" "$species" "$TOP_N" "$WINDOW" \
    > /dev/null

  echo "[run_confirm_positional_clade] ${species}: submitted"
done < "$ASSEMBLY_LIST"

echo "[run_confirm_positional_clade] done submitting."
