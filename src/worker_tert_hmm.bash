#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# worker_tert_hmm.bash <species list>
#
# The sensitive TERT search, for species where miniprot (src/worker_tert.bash)
# found nothing. One species (line $LSB_JOBINDEX + 1 of the list):
#   1. esl-translate: every ORF of >= 40 codons, both strands
#   2. hmmsearch of the TERT profile (src/build_tert_hmm.bash) and of
#      Telomerase_RBD (Pfam PF12009) against the ORFs
#      -> outputs/tert_hmm/<species>.{tert,rbd}.domtbl
# An ORF is one exon or less, so a gene's hits are split across ORFs;
# src/summarise_tert_hmm.py groups them by locus.
#
# Settings from the environment (exported by run_tert_hmm.bash):
# ESL_TRANSLATE, HMMSEARCH, TERT_HMM, RBD_HMM, OUT_DIR, THREADS.
# ----------------------------------------------------------------------

list="$1"
line=$(sed -n "$(( LSB_JOBINDEX + 1 ))p" "$list")
species=$(cut -f1 <<< "$line")
fasta=$(cut -f2 <<< "$line")
[[ -n "$species" && -n "$fasta" ]] || { echo "no species at index ${LSB_JOBINDEX} of ${list}" >&2; exit 1; }
echo "[worker_tert_hmm] ${species} ${fasta}"

tmp=$(mktemp -d "${TMPDIR:-/tmp}/tert_hmm.XXXXXX")
trap 'rm -rf "${tmp}"' EXIT

gzip -dcf "${fasta}" > "${tmp}/genome.fa"
"${ESL_TRANSLATE}" -l 40 "${tmp}/genome.fa" > "${tmp}/orfs.faa"
rm "${tmp}/genome.fa"

for model in tert rbd; do
    hmm="${TERT_HMM}"
    [[ "${model}" == "rbd" ]] && hmm="${RBD_HMM}"
    "${HMMSEARCH}" --cpu "${THREADS}" -E 1e-3 --domtblout "${tmp}/${model}.domtbl" \
        "${hmm}" "${tmp}/orfs.faa" > /dev/null
done
mv "${tmp}/tert.domtbl" "${OUT_DIR}/${species}.tert.domtbl"
mv "${tmp}/rbd.domtbl" "${OUT_DIR}/${species}.rbd.domtbl"
