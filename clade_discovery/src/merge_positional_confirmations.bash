#!/usr/bin/env bash
set -euo pipefail

# ./merge_positional_confirmations.bash
#
# Concatenates every clade_discovery/outputs/positional/summaries/<species>.tsv
# (one per species, written independently by parallel LSF jobs to avoid
# concurrent-write corruption on a shared file) into one combined table.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLADE_ROOT="$(dirname "$SCRIPT_DIR")"
SUMMARIES_DIR="${CLADE_ROOT}/outputs/positional/summaries"
OUT="${CLADE_ROOT}/outputs/positional_confirmations.tsv"

first=1
: > "$OUT"
for f in "$SUMMARIES_DIR"/*.tsv; do
  [[ -e "$f" ]] || continue
  if [[ "$first" == "1" ]]; then
    cat "$f" >> "$OUT"
    first=0
  else
    tail -n +2 "$f" >> "$OUT"
  fi
done

echo "[merge_positional_confirmations] $(wc -l < "$OUT") lines (incl. header) -> $OUT"
