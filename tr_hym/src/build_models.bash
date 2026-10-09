#!/usr/bin/env bash
set -euo pipefail
# ----------------------------------------------------------------------
# build_models.bash
#
# DNA profile HMMs of the published Hymenoptera TRs (Fajkus et al. 2023,
# NAR, Table S1 accessions; ../inputs/fajkus2023/), one per family with
# >= 3 TRs plus one of the smaller families pooled (models/<group>.fa):
# MAFFT alignment, then hmmbuild --dna. Concatenated into
# models/hymenoptera_tr.hmm for src/worker_search.bash.
# ----------------------------------------------------------------------
cd "$(dirname "$0")/../models"
MAFFT="/software/team301/mafft-7.525-with-extensions/scripts/mafft"
export MAFFT_BINARIES=/software/team301/mafft-7.525-with-extensions/binaries
HMMBUILD=/software/team301/hmmer-3.4/src/hmmbuild
: > hymenoptera_tr.hmm
for fa in *.fa; do
    [[ "$fa" == all_hymenoptera_tr.fa ]] && continue
    g="${fa%.fa}"
    "$MAFFT" --auto --quiet --thread 4 "$fa" > "${g}.aln.fa"
    "$HMMBUILD" --dna -n "TR_${g}" "${g}.hmm" "${g}.aln.fa" > /dev/null
    cat "${g}.hmm" >> hymenoptera_tr.hmm
    echo "${g}: $(grep -c '>' "$fa") TRs, $(grep -m1 '^LENG' "${g}.hmm" | awk '{print $2}') match states"
done
