#!/usr/bin/env bash
set -euo pipefail

# ./build_clade_assembly_list.bash <clade_name> <genus_list_file> [darwin_group1,darwin_group2,...]
#
# Builds a species<TAB>path assembly list for a target clade, scoped to
# an explicit genus list (produced by filtering classify_genera_ncbi.py's
# output TSV for the relevant NCBI group_name(s) - e.g. "spiders" for
# Araneae, "centipedes"+"millipedes" for Myriapoda, "beetles" for
# Coleoptera). Mirrors src/update_assembly_list.bash's organelle/cobiont
# SUSPECT screening (same naming-pattern + size-floor + plastid/mito
# keyword checks, same MIN_BYTES default) so this new clade's list starts
# from the same proven safety net rather than a naive re-implementation -
# see notes/organelle_and_cobiont_assembly_pitfall.md in the parent repo
# for why this matters (a real DToL gotcha, not hypothetical).
#
# Output: clade_discovery/inputs/<clade_name>_paths.txt

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLADE_ROOT="$(dirname "$SCRIPT_DIR")"
PROJECT_ROOT="$(dirname "$CLADE_ROOT")"

CLADE_NAME="${1:?usage: build_clade_assembly_list.bash <clade_name> <genus_list_file> [darwin_groups]}"
GENUS_LIST="${2:?usage: build_clade_assembly_list.bash <clade_name> <genus_list_file> [darwin_groups]}"
DARWIN_GROUPS="${3:-arthropods,insects}"

LATEST_TOOL="${LATEST_TOOL:-/software/team301/user/mb39/scripts/any_latest_assemblies.bash}"
MIN_BYTES="${MIN_BYTES:-5000000}"

OUT="${CLADE_ROOT}/inputs/${CLADE_NAME}_paths.txt"
SUSPECT="${CLADE_ROOT}/inputs/${CLADE_NAME}_suspect.txt"

echo "[build_clade_assembly_list] clade: $CLADE_NAME"
echo "[build_clade_assembly_list] genus list: $GENUS_LIST ($(wc -l < "$GENUS_LIST") genera)"
echo "[build_clade_assembly_list] darwin groups: $DARWIN_GROUPS"

IFS=',' read -ra DARWIN_GROUP_ARR <<< "$DARWIN_GROUPS" || true
RAW=$(mktemp)
bash "$LATEST_TOOL" "${DARWIN_GROUP_ARR[@]}" | sort > "$RAW"
echo "[build_clade_assembly_list] $(wc -l < "$RAW") species total across requested groups"

# Filter to species whose genus (first underscore-delimited token) is in
# the target genus list.
FILTERED=$(mktemp)
awk -F'\t' 'NR==FNR{genus[$1]=1; next} {
  split($1, a, "_"); if (a[1] in genus) print
}' "$GENUS_LIST" "$RAW" > "$FILTERED"
echo "[build_clade_assembly_list] $(wc -l < "$FILTERED") species matched target genera"

: > "$OUT"
: > "$SUSPECT"
while IFS=$'\t' read -r species path; do
  [[ -z "$species" ]] && continue
  reason=""
  curdir=$(basename "$(dirname "$path")")
  if echo "$path" | grep -qiE 'plastid|chloroplast|mito(chondri)?'; then
    reason="path mentions plastid/chloroplast/mito"
  elif ! echo "$curdir" | grep -qE '^[A-Za-z]+[0-9]+(\.hap[0-9]+)?\.[0-9]+$'; then
    reason="curated-directory name '${curdir}' doesn't match the plain assembly-version pattern"
  elif [[ -f "$path" ]]; then
    sz=$(stat -c%s "$path" 2>/dev/null || echo 0)
    if [[ "$sz" -gt 0 && "$sz" -lt "$MIN_BYTES" ]]; then
      reason="only ${sz} bytes (< ${MIN_BYTES})"
    fi
  fi
  if [[ -n "$reason" ]]; then
    echo -e "${species}\t${path}\t${reason}" >> "$SUSPECT"
  else
    echo -e "${species}\t${path}" >> "$OUT"
  fi
done < "$FILTERED"
sort -o "$OUT" "$OUT"

n_suspect=$(wc -l < "$SUSPECT")
echo "[build_clade_assembly_list] $(wc -l < "$OUT") species written to $OUT"
if [[ "$n_suspect" -gt 0 ]]; then
  echo "[build_clade_assembly_list] $n_suspect SUSPECT (excluded) - see $SUSPECT"
fi

rm -f "$RAW" "$FILTERED"
