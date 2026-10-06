#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# run_tidk_ends_candidates.bash
#
# After run_tidk_v03.bash: run tidk ends for each species' top
# N_CANDIDATES explore units on LSF (see worker_tidk_ends_candidates.bash),
# so the telomere can be called as the candidate found at the most
# chromosome ends rather than simply the most abundant explore unit.
#
# Reads species<TAB>path from: inputs/dtol_plant_paths.txt (or ASSEMBLY_LIST)
# Outputs: outputs/tidk_ends/candidates/<species>.<rank>.<unit>.ends.*
# Skips species that already have candidate output (unless FORCE=1), or
# that have no explore output yet.
# ----------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

ASSEMBLY_LIST="${ASSEMBLY_LIST:-${PROJECT_ROOT}/inputs/dtol_plant_paths.txt}"
EXPLORE_DIR="${EXPLORE_DIR:-${OUTPUT_DIR}/tidk_v0.3}"
CAND_DIR="${CAND_DIR:-${OUTPUT_DIR}/tidk_ends/candidates}"
LOGDIR="${OUTPUT_DIR}/logs/tidk_ends_candidates"
mkdir -p "$CAND_DIR" "$LOGDIR"

export TIDK="${TIDK:-${PROJECT_ROOT}/bin/tidk-0.3.0}"
export MIN_SEQ="${MIN_SEQ:-1000000}"
export N_CANDIDATES="${N_CANDIDATES:-5}"

memGB="${memGB:-4}"
QUEUE="${QUEUE:-normal}"
memMB=$(( memGB * 1024 ))

"${TIDK}" --version | grep -q "tidk 0.3" || { echo "[run_tidk_ends_candidates] ${TIDK} is not tidk 0.3.x" >&2; exit 1; }

submitted=0
while IFS=$'\t' read -r species fasta; do
    [[ -z "$species" || -z "$fasta" ]] && continue
    [[ -s "${EXPLORE_DIR}/${species}.explore.tsv" ]] || continue
    if compgen -G "${CAND_DIR}/${species}.1.*.ends.json" > /dev/null && [[ "${FORCE:-0}" != "1" ]]; then
        continue
    fi
    bsub \
      -J "tidkcand_${species}" \
      -q "$QUEUE" \
      -n 1 \
      -R "select[mem>${memMB}] rusage[mem=${memMB}]" \
      -M "$memMB" \
      -o "${LOGDIR}/${species}.%J.out" \
      -e "${LOGDIR}/${species}.%J.err" \
      bash "${SCRIPT_DIR}/worker_tidk_ends_candidates.bash" "$species" "$fasta" "$EXPLORE_DIR" "$CAND_DIR" \
      > /dev/null
    submitted=$(( submitted + 1 ))
done < "$ASSEMBLY_LIST"

echo "[run_tidk_ends_candidates] Submitted ${submitted} jobs (top ${N_CANDIDATES} units each)."
