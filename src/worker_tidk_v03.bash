#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# worker_tidk_v03.bash <species> <fasta> <explore_dir> <ends_dir>
#
# One species, run by run_tidk_v03.bash on LSF (or by hand):
#   1. tidk explore -> <explore_dir>/<species>.explore.tsv
#      (stderr, including strand-bias warnings, -> <species>.explore.err)
#   2. tidk ends with the top explore unit as --string (or discovering the
#      unit itself if explore found nothing)
#      -> <ends_dir>/<species>.ends.{tsv,bed,json}
#
# Settings come from the environment (exported by run_tidk_v03.bash):
# TIDK, MINLEN, MAXLEN, THRESH, DIST, MIN_SEQ.
# ----------------------------------------------------------------------

species="$1"
fasta="$2"
explore_dir="$3"
ends_dir="$4"

explore_out="${explore_dir}/${species}.explore.tsv"
explore_err="${explore_dir}/${species}.explore.err"

# write to a temporary file first, so a killed job doesn't leave a
# partial output that later runs would skip as done
"${TIDK}" explore \
    -m "${MINLEN}" -x "${MAXLEN}" \
    --threshold "${THRESH}" \
    --distance "${DIST}" \
    "${fasta}" > "${explore_out}.tmp" 2> "${explore_err}"
mv "${explore_out}.tmp" "${explore_out}"

unit=$(awk -F'\t' 'NR == 2 { print $1 }' "${explore_out}")
string_arg=()
if [[ -n "${unit}" ]]; then
    string_arg=(--string "${unit}")
fi

"${TIDK}" ends \
    "${string_arg[@]}" \
    --min-sequence-length "${MIN_SEQ}" \
    -o "${species}" -d "${ends_dir}" \
    "${fasta}" 2> "${ends_dir}/${species}.ends.err"
