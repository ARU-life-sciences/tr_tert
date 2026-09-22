#!/usr/bin/env bash
set -euo pipefail

# ./reciprocal_filter_candidates.bash <clade_name> <query_species> <subject_species> [evalue] [outdir]
#
# Step 4 of the clade_discovery method: blasts one species' candidate
# TR-like loci (predict_tr_candidate_loci.py's output - a few hundred
# to ~100k short context-extended sequences) against ANOTHER species'
# FULL GENOME. Most candidates are noise (near-chance short-tandem
# matches - see predict_tr_candidate_loci.py's docstring for why the
# raw candidate count is large) and their flanking context won't be
# conserved in an unrelated genome by chance; a genuine, functionally
# conserved locus (a real TR gene) should show a strong hit including
# flanking sequence, not just the repeat itself. This mirrors the
# Bombus paper's "mutual blastn searches (e-value threshold <= 1e-20)"
# step - same e-value default here.
#
# Submitted via LSF (this is a real compute job - candidate counts have
# been in the tens of thousands per species in practice, against a
# multi-hundred-Mb genome).

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLADE_ROOT="$(dirname "$SCRIPT_DIR")"

CLADE_NAME="${1:?usage: reciprocal_filter_candidates.bash <clade_name> <query_species> <subject_species> [evalue] [outdir]}"
QUERY_SP="${2:?usage: reciprocal_filter_candidates.bash <clade_name> <query_species> <subject_species> [evalue] [outdir]}"
SUBJECT_SP="${3:?usage: reciprocal_filter_candidates.bash <clade_name> <query_species> <subject_species> [evalue] [outdir]}"
EVALUE="${4:-1e-20}"
OUTDIR="${5:-${CLADE_ROOT}/outputs/reciprocal_filter}"
mkdir -p "$OUTDIR"

BLASTBIN="${BLASTBIN:-/software/team301/ncbi-blast-2.16.0+/bin}"
CANDIDATE_FA="${CLADE_ROOT}/outputs/candidate_loci/${QUERY_SP}_candidates.fa"
[[ -s "$CANDIDATE_FA" ]] || { echo "No candidates at $CANDIDATE_FA - run predict_tr_candidate_loci.py first" >&2; exit 1; }

SUBJECT_PATH=$(awk -F'\t' -v s="$SUBJECT_SP" '$1==s{print $2}' "${CLADE_ROOT}/inputs/${CLADE_NAME}_paths.txt")
[[ -n "$SUBJECT_PATH" ]] || { echo "Subject species $SUBJECT_SP not found in ${CLADE_NAME}_paths.txt" >&2; exit 1; }

QUEUE="${QUEUE:-normal}"
memGB="${memGB:-8}"
nCPU="${nCPU:-8}"
memMB=$(( memGB * 1024 ))

LOGDIR="${CLADE_ROOT}/outputs/logs/reciprocal_filter"
mkdir -p "$LOGDIR"

jobname="recip_${QUERY_SP}_vs_${SUBJECT_SP}"
out_hits="${OUTDIR}/${QUERY_SP}_vs_${SUBJECT_SP}.tsv"

echo "[reciprocal_filter_candidates] ${QUERY_SP} candidates vs ${SUBJECT_SP} genome (e<=${EVALUE})"

bsub \
  -J "$jobname" \
  -q "$QUEUE" \
  -n "$nCPU" \
  -R "span[hosts=1] select[mem>${memMB}] rusage[mem=${memMB}]" \
  -M "$memMB" \
  -o "${LOGDIR}/${jobname}.%J.out" \
  -e "${LOGDIR}/${jobname}.%J.err" \
  bash -lc "
    tmpdir=\$(mktemp -d -p \${TMPDIR:-/tmp} recip_${QUERY_SP}_${SUBJECT_SP}.XXXXXX)
    trap 'rm -rf \"\$tmpdir\"' EXIT
    fa=\"\$tmpdir/subject.fa\"
    if [[ '$SUBJECT_PATH' =~ \.gz\$ ]]; then gunzip -c '$SUBJECT_PATH' > \"\$fa\"; else cp '$SUBJECT_PATH' \"\$fa\"; fi
    '$BLASTBIN/makeblastdb' -in \"\$fa\" -dbtype nucl -out \"\$tmpdir/db\"
    '$BLASTBIN/blastn' -query '$CANDIDATE_FA' -db \"\$tmpdir/db\" \
      -num_threads $nCPU -evalue $EVALUE \
      -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
      -out '${out_hits}.tmp'
    mv '${out_hits}.tmp' '$out_hits'
  "

echo "[reciprocal_filter_candidates] submitted, output -> $out_hits"
