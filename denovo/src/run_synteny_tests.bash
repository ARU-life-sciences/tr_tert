#!/usr/bin/env bash
set -euo pipefail
# Runs synteny_test.py on a pairs file (distance<TAB>A<TAB>B; both TRs known
# from a clade's template_check.tsv). Defaults: synteny_pairs.tsv and
# ../tr_hym -> outputs/synteny_tests.tsv; override with PAIRS, CLADE, OUT.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CLADE="${CLADE:-${ROOT}/../tr_hym}"
TC="${CLADE}/outputs/template_check.tsv"
SP="${CLADE}/species.tsv"
PAIRS="${PAIRS:-${ROOT}/synteny_pairs.tsv}"
out="${OUT:-${ROOT}/outputs/synteny_tests.tsv}"
echo -e "distance\tA\tB\tanchors\tB_seq\tinterval\tinterval_kb\tB_TR_inside\twindows\tB_TR_window" > "$out"
tail -n +2 "$PAIRS" | while IFS=$'\t' read -r dist A B; do
    gA=$(awk -F'\t' -v s="$A" '$1==s{print $2}' "$SP"); gB=$(awk -F'\t' -v s="$B" '$1==s{print $2}' "$SP")
    lA=$(awk -F'\t' -v s="$A" '$1==s{print $5}' "$TC"); lB=$(awk -F'\t' -v s="$B" '$1==s{print $5}' "$TC")
    uB=$(awk -F'\t' -v s="$B" '$1==s{print $3}' "$TC")
    res=$(python3 "${ROOT}/src/synteny_test.py" "$gA" "$lA" "$gB" "$uB" "$lB" "${ROOT}/outputs/synteny/${A}__${B}")
    echo -e "${dist}\t${A}\t${B}\t${res}" >> "$out"
done
