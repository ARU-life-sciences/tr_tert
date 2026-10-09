#!/usr/bin/env bash
set -euo pipefail
# ----------------------------------------------------------------------
# worker_search.bash <species list>
#
# One species (line $LSB_JOBINDEX + 1 of the list: species<TAB>fasta...):
# nhmmer of the Hymenoptera TR profile HMMs (models/hymenoptera_tr.hmm)
# (--dna: a T2T assembly can start with kilobases of G-free telomere,
# which defeats nhmmer's alphabet guess)
# against the assembly, and the hits (E <= 1e-5) with 300 bp either side,
# in the TR's orientation:
#   outputs/<species>.nhmmer.tbl, outputs/<species>.hits.fa
# ----------------------------------------------------------------------
list="$1"
line=$(sed -n "$(( LSB_JOBINDEX + 1 ))p" "$list")
species=$(cut -f1 <<< "$line"); fasta=$(cut -f2 <<< "$line")
echo "[worker_search] ${species} ${fasta}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
H=/software/team301/hmmer-3.4
tmp=$(mktemp -d "${TMPDIR:-/tmp}/trhym.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
gzip -dcf "$fasta" > "$tmp/g.fa"
"$H/easel/miniapps/esl-sfetch" --index "$tmp/g.fa" > /dev/null
"$H/src/nhmmer" --dna --cpu 8 -E 1e-5 --tblout "$tmp/hits.tbl" "$ROOT/models/hymenoptera_tr.hmm" "$tmp/g.fa" > /dev/null
# seq, from, to (from > to on the minus strand), model, E-value; no hits
# (grep exits 1) is a result, not a failure
{ grep -v '^#' "$tmp/hits.tbl" || true; } | awk '{print $1, $7, $8, $3, $13}' | while read -r seq a b model e; do
    if (( a < b )); then s=$(( a > 300 ? a - 300 : 1 )); t=$(( b + 300 )); else s=$(( a + 300 )); t=$(( b > 300 ? b - 300 : 1 )); fi
    "$H/easel/miniapps/esl-sfetch" -n "${seq}:${a}-${b}|${model}|${e}" -c "${s}..${t}" "$tmp/g.fa" "$seq" 2>/dev/null || true
done > "$tmp/hits.fa"
mv "$tmp/hits.tbl" "$ROOT/outputs/${species}.nhmmer.tbl"
mv "$tmp/hits.fa" "$ROOT/outputs/${species}.hits.fa"
