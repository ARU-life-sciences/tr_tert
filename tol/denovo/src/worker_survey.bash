#!/usr/bin/env bash
set -euo pipefail
# ----------------------------------------------------------------------
# worker_survey.bash <survey_pairs.tsv>
#
# One TR-synteny survey pair (line $LSB_JOBINDEX + 1): is B's known TR in
# the interval projected from A's TR, by nucleotide anchors
# (synteny_test.py) and by BUSCO-gene anchors (busco_synteny.py)?
#   -> outputs/survey/<n>.tsv: clade level A B nuc_anchors nuc_inside
#                              busco_anchors busco_inside
# ----------------------------------------------------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
i=${LSB_JOBINDEX}
IFS=$'\t' read -r clade level A lA fA B lB fB uB lin <<< "$(sed -n "$(( i + 1 ))p" "$1")"
d="${ROOT}/outputs/survey/work/${A}__${B}"
mkdir -p "$d"
nuc=$(python3 "${ROOT}/src/synteny_test.py" "$fA" "$lA" "$fB" "$uB" "$lB" "$d/nuc" || echo "error")
busco=$(python3 "${ROOT}/src/busco_synteny.py" "$fA" "$lA" "$fB" "$uB" "$lB" "$d/busco" - - "$lin" || echo "error")
na=$(cut -f1 <<< "$nuc"); ni=$(cut -f5 <<< "$nuc")
ba=$(cut -f1 <<< "$busco"); bi=$(cut -f5 <<< "$busco")
echo -e "${clade}\t${level}\t${A}\t${B}\t${na}\t${ni}\t${ba}\t${bi}" > "${ROOT}/outputs/survey/${i}.tsv"
rm -rf "$d"
