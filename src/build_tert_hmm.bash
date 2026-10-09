#!/usr/bin/env bash
set -euo pipefail
# ----------------------------------------------------------------------
# build_tert_hmm.bash <proteins.faa> <out prefix>
#
# A TERT profile HMM from the reference panel plus the DToL TERTs found by
# miniprot (inputs/tert_panel_round3.faa): MAFFT alignment, then hmmbuild.
# For the sensitive search of species where pairwise alignment found
# nothing (src/worker_tert_hmm.bash).
# ----------------------------------------------------------------------
in="$1"
out="$2"
MAFFT_BINARIES=/software/team301/mafft-7.525-with-extensions/binaries \
  /software/team301/mafft-7.525-with-extensions/scripts/mafft --auto --thread 8 --quiet "$in" > "${out}.aln.faa"
/software/team301/hmmer-3.4/src/hmmbuild --amino -n TERT_panel "${out}.hmm" "${out}.aln.faa" > "${out}.hmmbuild.log"
grep -m1 "^LENG" "${out}.hmm"
