#!/usr/bin/env bash
set -euo pipefail

# ----------------------------------------------------------------------
# worker_tert.bash <species list>
#
# One species (line $LSB_JOBINDEX + 1 of the list, after its header):
#   1. miniprot: spliced alignment of the TERT reference panel
#      (inputs/tert_panel.faa) to the assembly, with translations
#      -> outputs/tert/<species>.miniprot.paf
#   2. hmmsearch of the Telomerase_RBD domain (Pfam PF12009) against the
#      translated alignments, which tells TERT from the reverse
#      transcriptases of retrotransposons (which share its RT domain)
#      -> outputs/tert/<species>.rbd.domtbl
# src/summarise_tert.py turns these into one call per species.
#
# Settings from the environment (exported by run_tert.bash): MINIPROT,
# HMMSEARCH, PANEL, RBD_HMM, OUT_DIR, THREADS.
# ----------------------------------------------------------------------

list="$1"
line=$(sed -n "$(( LSB_JOBINDEX + 1 ))p" "$list")
species=$(cut -f1 <<< "$line")
fasta=$(cut -f2 <<< "$line")
[[ -n "$species" && -n "$fasta" ]] || { echo "no species at index ${LSB_JOBINDEX} of ${list}" >&2; exit 1; }
echo "[worker_tert] ${species} ${fasta}"

paf="${OUT_DIR}/${species}.miniprot.paf"
domtbl="${OUT_DIR}/${species}.rbd.domtbl"

tmp=$(mktemp -d "${TMPDIR:-/tmp}/tert.XXXXXX")
trap 'rm -rf "${tmp}"' EXIT

# -I: intron size from genome size; --outn: up to 20 alignments per query,
# so paralogues and retrotransposon RT hits are seen and can be told apart
"${MINIPROT}" -t "${THREADS}" -I --outn=20 --trans "${fasta}" "${PANEL}" > "${tmp}/out.paf"

# translations are "##STA <protein>" lines after each alignment
awk -F'\t' '
    /^##STA/ { if (name != "") { print ">" name; print $2 } ; name = ""; next }
    /^#/ { next }
    { name = $6 ":" $8 "-" $9 ":" $5 "|" $1 }
' "${tmp}/out.paf" > "${tmp}/trans.faa"

if [[ -s "${tmp}/trans.faa" ]]; then
    "${HMMSEARCH}" --cpu "${THREADS}" --domtblout "${tmp}/rbd.domtbl" -E 1e-3 \
        "${RBD_HMM}" "${tmp}/trans.faa" > /dev/null
else
    # no TERT-like alignment at all (e.g. Diptera): an empty result
    echo "# no miniprot alignments" > "${tmp}/rbd.domtbl"
fi

mv "${tmp}/out.paf" "${paf}"
mv "${tmp}/rbd.domtbl" "${domtbl}"
