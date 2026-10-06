#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# run_identical_set_array_structure_bsub.bash submit|merge <variant_positions.tsv> <outdir> [window_bp]
#
# LSF counterpart to tr_identical_set_array_structure.py, which walks
# every main chromosome end of each species in turn and is slow when run
# over many genomes serially (each genome is decompressed and indexed).
#
#   submit: split the variant-positions TSV (from
#           find_identical_set_second_candidates.py) by species and submit
#           one walk job per species -> <outdir>/per_species/<species>.tsv
#   merge:  concatenate the per-species results (one header) into
#           <outdir>/identical_set_array_structure.tsv, and list any
#           species without a result
#
# Skips species that already have a result (unless FORCE=1).
# ----------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mode="$1"
vp="$2"
outdir="$3"
window="${4:-5000}"

per_species="${outdir}/per_species"
inputs="${outdir}/inputs"
logdir="${outdir}/logs"
mkdir -p "$per_species" "$inputs" "$logdir"

memGB="${memGB:-8}"
memMB=$(( memGB * 1024 ))

species_list() {
    awk -F'\t' 'NR > 1 { print $1 }' "$vp" | sort -u
}

case "$mode" in
submit)
    submitted=0
    for species in $(species_list); do
        out="${per_species}/${species}.tsv"
        if [[ -s "$out" && "${FORCE:-0}" != "1" ]]; then
            continue
        fi
        # this species' rows, with the header
        awk -F'\t' -v sp="$species" 'NR == 1 || $1 == sp' "$vp" > "${inputs}/${species}.tsv"
        bsub -J "idwalk_${species}" -n 1 -q normal \
          -R "select[mem>${memMB}] rusage[mem=${memMB}]" \
          -M "${memMB}" \
          -o "${logdir}/${species}.%J.out" \
          -e "${logdir}/${species}.%J.err" \
          "python3 '${SCRIPT_DIR}/tr_identical_set_array_structure.py' '${inputs}/${species}.tsv' '${window}' '${out}'" \
          > /dev/null
        submitted=$(( submitted + 1 ))
    done
    echo "[run_identical_set_array_structure_bsub] Submitted ${submitted} jobs. Outputs: ${per_species}"
    ;;
merge)
    merged="${outdir}/identical_set_array_structure.tsv"
    first=1
    : > "$merged"
    for species in $(species_list); do
        out="${per_species}/${species}.tsv"
        if [[ ! -s "$out" ]]; then
            echo "[run_identical_set_array_structure_bsub] no result yet: ${species}" >&2
            continue
        fi
        if (( first )); then
            cat "$out" >> "$merged"
            first=0
        else
            tail -n +2 "$out" >> "$merged"
        fi
    done
    echo "[run_identical_set_array_structure_bsub] Merged -> ${merged}"
    ;;
*)
    echo "usage: $0 submit|merge <variant_positions.tsv> <outdir> [window_bp]" >&2
    exit 1
    ;;
esac
