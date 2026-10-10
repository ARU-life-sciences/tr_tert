#!/usr/bin/env bash
set -euo pipefail
# ----------------------------------------------------------------------
# worker_tol.bash <species list>
#
# ToL-wide TR search, one genome (line $LSB_JOBINDEX + 1 of the list):
# every TR model we have, whatever clade it came from -
#   models/library.hmm  nhmmer: Hymenoptera family/lineage HMMs (incl.
#                       DToL seeds), Pimplinae, Lepidoptera family HMMs,
#                       the tr_tert plant TR HMM
#   models/library.cm   cmsearch: Rfam vertebrate, ciliate and yeast TR CMs
# Hits (E <= 1e-5) with 300 bp flanks in TR orientation ->
# outputs/<species>.{nhmmer.tbl,hits.fa} (hits.fa headers seq:a-b|model|E),
# for template_check.py (copy-aware: PICK_TEMPLATE=1).
# ----------------------------------------------------------------------
# CLADE_DIR: the round directory (models/library.{hmm,cm}, outputs/); default tr_tol
ROOT="${CLADE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
line=$(sed -n "$(( LSB_JOBINDEX + 1 ))p" "$1")
species=$(cut -f1 <<< "$line"); fasta=$(cut -f2 <<< "$line")
echo "[worker_tol] ${species} ${fasta}"
H=/software/team301/hmmer-3.4
I=/software/team301/infernal-1.1.4-linux-intel-gcc/binaries
tmp=$(mktemp -d "${TMPDIR:-/tmp}/trtol.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
# NUL bytes occur in some curated assemblies (fSebMen2.1, mMyoGla1.1)
gzip -dcf "$fasta" | tr -d '\000' > "$tmp/g.fa"
"$H/easel/miniapps/esl-sfetch" --index "$tmp/g.fa" > /dev/null
# --dna: T2T assemblies can start with kilobases of G-free telomere
"$H/src/nhmmer" --dna --cpu 8 -E 1e-5 --tblout "$tmp/n.tbl" "$ROOT/models/library.hmm" "$tmp/g.fa" > /dev/null
# CMs only if the round has them (walking rounds search seed HMMs only)
if [[ -s "$ROOT/models/library.cm" ]]; then
    "$I/cmsearch" --cpu 8 -E 1e-5 --tblout "$tmp/cm.tbl" "$ROOT/models/library.cm" "$tmp/g.fa" > /dev/null
else
    : > "$tmp/cm.tbl"
fi
# common columns: target, -, model, -, -, -, from, to, ..., E (13)
{ grep -v '^#' "$tmp/n.tbl" || true
  awk '!/^#/ {print $1, "-", $3, "-", "-", "-", $8, $9, "-", "-", "-", "-", $16}' "$tmp/cm.tbl"
} > "$tmp/hits.tbl"
awk '{print $1, $7, $8, $3, $13}' "$tmp/hits.tbl" | while read -r seq a b model e; do
    if (( a < b )); then s=$(( a > 300 ? a - 300 : 1 )); t=$(( b + 300 )); else s=$(( a + 300 )); t=$(( b > 300 ? b - 300 : 1 )); fi
    "$H/easel/miniapps/esl-sfetch" -n "${seq}:${a}-${b}|${model}|${e}" -c "${s}..${t}" "$tmp/g.fa" "$seq" 2>/dev/null || true
done > "$tmp/hits.fa"
mv "$tmp/hits.tbl" "$ROOT/outputs/${species}.nhmmer.tbl"
mv "$tmp/hits.fa" "$ROOT/outputs/${species}.hits.fa"
