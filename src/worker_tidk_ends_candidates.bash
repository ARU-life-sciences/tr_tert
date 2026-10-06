#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# worker_tidk_ends_candidates.bash <species> <fasta> <explore_dir> <cand_dir>
#
# Run tidk ends for each of the top N_CANDIDATES units from a species'
# tidk 0.3.0 explore output (<explore_dir>/<species>.explore.tsv), writing
# <cand_dir>/<species>.<rank>.<unit>.ends.{tsv,bed,json}.
#
# The top explore unit isn't always the telomere: a long subterminal
# satellite can have more copies near the chromosome ends. Calling ends
# for several candidates lets the telomere be picked as the unit actually
# found at the most chromosome ends (src/tidk_v03_summary.py).
#
# Settings from the environment (exported by run_tidk_ends_candidates.bash):
# TIDK, MIN_SEQ, N_CANDIDATES.
# ----------------------------------------------------------------------

species="$1"
fasta="$2"
explore_dir="$3"
cand_dir="$4"

rank=0
awk -F'\t' -v n="${N_CANDIDATES}" 'NR > 1 && NR <= n + 1 { print $1 }' \
    "${explore_dir}/${species}.explore.tsv" |
while read -r unit; do
    rank=$(( rank + 1 ))
    "${TIDK}" ends \
        --string "${unit}" \
        --min-sequence-length "${MIN_SEQ}" \
        -o "${species}.${rank}.${unit}" -d "${cand_dir}" \
        "${fasta}" 2> /dev/null
done
