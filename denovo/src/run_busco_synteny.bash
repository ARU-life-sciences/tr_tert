#!/usr/bin/env bash
set -euo pipefail
# busco_synteny.py over a pairs file: distance A A_TR B unit_B B_TR(or none) hmm cm(or -)
# genomes from ../tr_hym/species.tsv. Usage: run_busco_synteny.bash <pairs.tsv> <out.tsv>
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SP="${ROOT}/../tr_hym/species.tsv"
echo -e "distance\tA\tB\tanchors\tB_seq\tinterval\tinterval_kb\tB_TR_inside\tbest_hit\tbest_E\tbest_template\tn_hits_with_template" > "$2"
tail -n +2 "$1" | while IFS=$'\t' read -r dist A lA B uB lB hmm cm; do
    gA=$(awk -F'\t' -v s="$A" '$1==s{print $2}' "$SP"); gB=$(awk -F'\t' -v s="$B" '$1==s{print $2}' "$SP")
    res=$(python3 "${ROOT}/src/busco_synteny.py" "$gA" "$lA" "$gB" "$uB" "$lB" "${ROOT}/outputs/busco_synteny/${A}__${B}" "$hmm" "$cm" 2>"${ROOT}/outputs/busco_synteny/err.${A}__${B}" || echo "error")
    echo -e "${dist}\t${A}\t${B}\t${res}" >> "$2"
done
