#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# run_tidk_v03.bash
#
# Run tidk 0.3.0 "explore" and "ends" across assemblies using LSF, one
# job per species (see worker_tidk_v03.bash). Taxon-agnostic: any
# species<TAB>fasta list works, plants are just the default.
#
# Reads species<TAB>path from: inputs/dtol_plant_paths.txt (or ASSEMBLY_LIST)
# Outputs: outputs/tidk_v0.3/<species>.explore.tsv
#          outputs/tidk_ends/<species>.ends.{tsv,bed,json}
# Skips species with an existing ends JSON (unless FORCE=1).
#
# The v0.2.7 results in outputs/tidk/ are left untouched: they back the
# existing notes. 0.3.0 differs in ways that matter downstream: repeats
# are collapsed to their primitive unit (no multimers), counts are copies
# rather than run counts, the output has four columns
# (canonical_repeat_unit, copies, copies_as_unit, copies_as_revcomp),
# and runs tolerate sequencing errors by default. See the tidk CHANGELOG.
# ----------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_ROOT}/inputs"
OUTPUT_DIR="${PROJECT_ROOT}/outputs"

ASSEMBLY_LIST="${ASSEMBLY_LIST:-${INPUT_DIR}/dtol_plant_paths.txt}"

EXPLORE_DIR="${EXPLORE_DIR:-${OUTPUT_DIR}/tidk_v0.3}"
ENDS_DIR="${ENDS_DIR:-${OUTPUT_DIR}/tidk_ends}"
LOGDIR="${OUTPUT_DIR}/logs/tidk_v0.3"
mkdir -p "$EXPLORE_DIR" "$ENDS_DIR" "$LOGDIR"

# a fixed 0.3.0 binary, so results don't change if tidk is rebuilt
export TIDK="${TIDK:-${PROJECT_ROOT}/bin/tidk-0.3.0}"

# tidk explore parameters: the same as the v0.2.7 run, so the two compare
export MINLEN="${MINLEN:-4}"      # -m
export MAXLEN="${MAXLEN:-30}"     # -x
export THRESH="${THRESH:-100}"    # --threshold
export DIST="${DIST:-0.01}"       # --distance
# tidk ends: skip unplaced scaffolds and organelles
export MIN_SEQ="${MIN_SEQ:-1000000}"

nCPU="${nCPU:-8}"
memGB="${memGB:-8}"
QUEUE="${QUEUE:-normal}"
memMB=$(( memGB * 1024 ))

"${TIDK}" --version | grep -q "tidk 0.3" || { echo "[run_tidk_v03] ${TIDK} is not tidk 0.3.x" >&2; exit 1; }

echo "[run_tidk_v03] assembly list: $ASSEMBLY_LIST"
echo "[run_tidk_v03] outputs:       $EXPLORE_DIR, $ENDS_DIR"
echo "[run_tidk_v03] tidk:          $TIDK ($("$TIDK" --version))"
echo "[run_tidk_v03] explore:       -m $MINLEN -x $MAXLEN --threshold $THRESH --distance $DIST"
echo "[run_tidk_v03] ends:          --min-sequence-length $MIN_SEQ"
echo

submitted=0
while IFS=$'\t' read -r species fasta; do
    [[ -z "$species" || -z "$fasta" ]] && continue

    if [[ -s "${ENDS_DIR}/${species}.ends.json" && "${FORCE:-0}" != "1" ]]; then
        continue
    fi

    bsub \
      -J "tidk03_${species}" \
      -q "$QUEUE" \
      -n "$nCPU" \
      -R "span[hosts=1] select[mem>${memMB}] rusage[mem=${memMB}]" \
      -M "$memMB" \
      -o "${LOGDIR}/${species}.%J.out" \
      -e "${LOGDIR}/${species}.%J.err" \
      bash "${SCRIPT_DIR}/worker_tidk_v03.bash" "$species" "$fasta" "$EXPLORE_DIR" "$ENDS_DIR" \
      > /dev/null
    submitted=$(( submitted + 1 ))
done < "$ASSEMBLY_LIST"

echo "[run_tidk_v03] Submitted ${submitted} jobs."
