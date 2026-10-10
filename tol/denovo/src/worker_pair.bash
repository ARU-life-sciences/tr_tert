#!/usr/bin/env bash
set -euo pipefail
# ----------------------------------------------------------------------
# worker_pair.bash <pairs.tsv>
#
# One pair (line $LSB_JOBINDEX + 1 of pairs.tsv): src/covariation.py with
# A's candidates against B -> outputs/pairs/<A>__<B>/. If ../tr_hym found
# A's TR (template_check.tsv, E <= 1e-10), its locus is passed in so the
# summary says whether co-variation recovers it, and at what rank.
# ----------------------------------------------------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
line=$(sed -n "$(( LSB_JOBINDEX + 1 ))p" "$1")
IFS=$'\t' read -r genus family order A uA fA B uB fB <<< "$line"
out="${ROOT}/outputs/pairs/${A}__${B}"
mkdir -p "$out"
known=()
tc="${ROOT}/../tr_hym/outputs/template_check.tsv"
loc=$(awk -F'\t' -v s="$A" '$1==s && $5 ~ /:/ {print $5}' "$tc" 2>/dev/null || true)
if [[ -n "$loc" ]]; then
    echo "$loc" | awk -F'[:-]' '{print $1"\t"$2"\t"$3}' > "${out}/known_tr.tsv"
    known=("${out}/known_tr.tsv")
fi
echo "[worker_pair] ${genus}: ${A} (${uA}) vs ${B} (${uB})"
python3 "${ROOT}/src/covariation.py" "$fA" "$uA" "$fB" "$uB" "$out" "${known[@]}"
