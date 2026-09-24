#!/usr/bin/env bash
set -euo pipefail

# ./worker_reciprocal_filter.bash <candidate_fa> <subject_path> <evalue> <out_hits>
#
# Plain worker (see worker_snrna_search.bash for why this pattern is
# used instead of an inline bsub 'bash -lc "..."' string - the latter
# broke on nested quoting around `trap`).

candidate_fa="$1"
subject_path="$2"
evalue="$3"
out_hits="$4"

BLASTBIN="${BLASTBIN:-/software/team301/ncbi-blast-2.16.0+/bin}"
nCPU="${nCPU:-8}"

tmpdir=$(mktemp -d -p "${TMPDIR:-/tmp}" recip.XXXXXX)
trap 'rm -rf "$tmpdir"' EXIT

fa="${tmpdir}/subject.fa"
if [[ "$subject_path" =~ \.gz$ ]]; then
  gunzip -c "$subject_path" > "$fa"
else
  cp "$subject_path" "$fa"
fi

"$BLASTBIN/makeblastdb" -in "$fa" -dbtype nucl -out "${tmpdir}/db"
"$BLASTBIN/blastn" -query "$candidate_fa" -db "${tmpdir}/db" \
  -num_threads "$nCPU" -evalue "$evalue" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -out "${out_hits}.tmp"
mv "${out_hits}.tmp" "$out_hits"
