#!/usr/bin/env bash
set -euo pipefail

# ./run_snrna_search_clade.bash <clade_name> <species> [more species...]
#
# Runs cmsearch (Infernal) against Rfam's U1/U2/U3/U4/U5/U6/U6atac/
# RNase_MRP/Arthropod_7SK covariance models (clade_discovery/inputs/rfam/
# snRNA_promoter_models.cm - fetched from Rfam, see notes/plan.md) for
# one or more species. These are the RNAPIII type-3-promoter snRNA
# families used as anchors to characterise a species-specific USE/TATA
# promoter consensus (extract_snrna_promoters.py + build_promoter_consensus.py) -
# the unknown TR is presumed to share this promoter architecture, which
# is the whole premise of this clade_discovery effort.
#
# Submitted via LSF - one cmsearch run against a full genome across 9
# models took ~15 minutes wall-clock on 8 CPUs for a ~300Mb beetle
# genome in testing; budget accordingly for larger genomes.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLADE_ROOT="$(dirname "$SCRIPT_DIR")"

CLADE_NAME="${1:?usage: run_snrna_search_clade.bash <clade_name> <species> [more species...]}"
shift
SPECIES_LIST=("$@")
[[ ${#SPECIES_LIST[@]} -gt 0 ]] || { echo "No species given" >&2; exit 1; }

CM="${CLADE_ROOT}/inputs/rfam/snRNA_promoter_models.cm"
[[ -s "$CM" ]] || { echo "No CM file at $CM" >&2; exit 1; }
WORKER="${SCRIPT_DIR}/worker_snrna_search.bash"

ASSEMBLY_LIST="${CLADE_ROOT}/inputs/${CLADE_NAME}_paths.txt"

OUTDIR="${CLADE_ROOT}/outputs/snrna_search"
LOGDIR="${CLADE_ROOT}/outputs/logs/snrna_search"
mkdir -p "$OUTDIR" "$LOGDIR"

QUEUE="${QUEUE:-normal}"
memGB="${memGB:-8}"
nCPU="${nCPU:-8}"
memMB=$(( memGB * 1024 ))
FORCE="${FORCE:-0}"

for species in "${SPECIES_LIST[@]}"; do
  fasta=$(awk -F'\t' -v s="$species" '$1==s{print $2}' "$ASSEMBLY_LIST")
  if [[ -z "$fasta" ]]; then
    echo "[run_snrna_search_clade] ${species}: not found in ${CLADE_NAME}_paths.txt, skipping" >&2
    continue
  fi
  out_tbl="${OUTDIR}/${species}_snrna.tbl"
  if [[ "$FORCE" != "1" && -s "$out_tbl" ]]; then
    echo "[run_snrna_search_clade] ${species}: exists, skipping."
    continue
  fi

  jobname="snrna_${CLADE_NAME}_${species}"
  bsub \
    -J "$jobname" \
    -q "$QUEUE" \
    -n "$nCPU" \
    -R "span[hosts=1] select[mem>${memMB}] rusage[mem=${memMB}]" \
    -M "$memMB" \
    -o "${LOGDIR}/${species}.%J.out" \
    -e "${LOGDIR}/${species}.%J.err" \
    -env "all, nCPU=$nCPU" \
    bash "$WORKER" "$species" "$fasta" "$CM" "$out_tbl"
  echo "[run_snrna_search_clade] ${species}: submitted"
done
