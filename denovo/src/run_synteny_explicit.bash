#!/usr/bin/env bash
set -euo pipefail
# synteny_test.py over a pairs file with explicit loci and HMMs:
# distance A A_TR B unit_B B_TR(or none) hmm; genomes from ../tr_hym/species.tsv
# Usage: run_synteny_explicit.bash <pairs.tsv> <out.tsv>
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SP="${ROOT}/../tr_hym/species.tsv"
echo -e "distance\tA\tB\tanchors\tB_seq\tinterval\tinterval_kb\tB_TR_inside\twindows\tB_TR_window\thmm_best\thmm_best_template" > "$2"
tail -n +2 "$1" | while IFS=$'\t' read -r dist A lA B uB lB hmm; do
    gA=$(awk -F'\t' -v s="$A" '$1==s{print $2}' "$SP"); gB=$(awk -F'\t' -v s="$B" '$1==s{print $2}' "$SP")
    res=$(python3 "${ROOT}/src/synteny_test.py" "$gA" "$lA" "$gB" "$uB" "$lB" "${ROOT}/outputs/synteny/${A}__${B}" "$hmm" || echo "error")
    echo -e "${dist}\t${A}\t${B}\t${res}" >> "$2"
done
