#!/usr/bin/env bash
set -euo pipefail

# ./worker_snrna_search.bash <species> <fasta_path> <cm_path> <out_tbl>
#
# Plain worker script (no nested-quoting issues, unlike an earlier
# version of run_snrna_search_clade.bash that tried to inline this as a
# bsub 'bash -lc "..."' string and broke on nested quoting around
# `trap` - matches the parent pipeline's own worker_tblastn_tert.bash
# pattern instead).

species="$1"
fasta="$2"
cm="$3"
out_tbl="$4"

CMSEARCH="/software/team301/infernal-1.1.4-linux-intel-gcc/binaries/cmsearch"
nCPU="${nCPU:-8}"

tmpdir=$(mktemp -d -p "${TMPDIR:-/tmp}" "snrna_${species}.XXXXXX")
trap 'rm -rf "$tmpdir"' EXIT

fa="${tmpdir}/genome.fa"
if [[ "$fasta" =~ \.gz$ ]]; then
  gunzip -c "$fasta" > "$fa"
else
  cp "$fasta" "$fa"
fi

"$CMSEARCH" --cpu "$nCPU" --tblout "${out_tbl}.tmp" "$cm" "$fa"
mv "${out_tbl}.tmp" "$out_tbl"
