#!/usr/bin/env bash
set -euo pipefail

# ./run_tidk_clade.bash <clade_name>
#
# Stage 0 for a new clade: telomere motif discovery via `tidk explore`
# (de novo repeat discovery - no known motif required), across every
# species in clade_discovery/inputs/<clade_name>_paths.txt. Adapted from
# the parent pipeline's src/run_tidk.bash, kept as a separate copy
# (rather than parameterizing the original via env vars) so plant and
# clade-discovery runs can never collide on output paths, and so this
# component stays self-contained.
#
# This is the necessary prerequisite for everything downstream: motif
# VARIATION within a clade is what makes the reciprocal-filtering
# candidate-TR-locus search tractable at all (see
# clade_discovery/notes/plan.md) - without knowing what's actually out
# there, clade selection/prioritization is a guess.
#
# Outputs: clade_discovery/outputs/tidk/<species>.tidk.tsv
# Skips species with an existing output (unless FORCE=1).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLADE_ROOT="$(dirname "$SCRIPT_DIR")"

CLADE_NAME="${1:?usage: run_tidk_clade.bash <clade_name>}"
ASSEMBLY_LIST="${CLADE_ROOT}/inputs/${CLADE_NAME}_paths.txt"
[[ -s "$ASSEMBLY_LIST" ]] || { echo "No assembly list at $ASSEMBLY_LIST - run build_clade_assembly_list.bash first" >&2; exit 1; }

OUTDIR="${CLADE_ROOT}/outputs/tidk"
LOGDIR="${CLADE_ROOT}/outputs/logs/tidk"
mkdir -p "$OUTDIR" "$LOGDIR"

TIDK="${TIDK:-tidk}"
MINLEN="${MINLEN:-4}"      # -m
MAXLEN="${MAXLEN:-30}"     # -x
THRESH="${THRESH:-100}"    # --threshold (tidk default)
DIST="${DIST:-0.01}"       # --distance
FORCE="${FORCE:-0}"

QUEUE="${QUEUE:-normal}"  # 'small' queue no longer exists on this cluster
memGB="${memGB:-8}"
nCPU="${nCPU:-8}"
memMB=$(( memGB * 1024 ))

echo "[run_tidk_clade] clade: $CLADE_NAME"
echo "[run_tidk_clade] assembly list: $ASSEMBLY_LIST ($(wc -l < "$ASSEMBLY_LIST") species)"
echo "[run_tidk_clade] output: $OUTDIR"
echo "[run_tidk_clade] tidk params: -m $MINLEN -x $MAXLEN --threshold $THRESH --distance $DIST"

while IFS=$'\t' read -r species fasta; do
  [[ -z "$species" || -z "$fasta" ]] && continue
  out="${OUTDIR}/${species}.tidk.tsv"
  if [[ "$FORCE" != "1" && -s "$out" ]]; then
    echo "[run_tidk_clade] ${species}: exists, skipping."
    continue
  fi

  log_out="${LOGDIR}/${species}.%J.out"
  log_err="${LOGDIR}/${species}.%J.err"

  bsub \
    -J "tidk_${CLADE_NAME}_${species}" \
    -q "$QUEUE" \
    -n "$nCPU" \
    -R "span[hosts=1] select[mem>${memMB}] rusage[mem=${memMB}]" \
    -M "$memMB" \
    -o "$log_out" \
    -e "$log_err" \
    bash -lc "
      ${TIDK} explore \
          -m ${MINLEN} -x ${MAXLEN} \
          --threshold ${THRESH} \
          --distance ${DIST} \
          --log \
          '${fasta}' \
          > '${out}'
    "

  echo "[run_tidk_clade] ${species}: submitted"
done < "$ASSEMBLY_LIST"

echo "[run_tidk_clade] done submitting."
