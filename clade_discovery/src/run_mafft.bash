#!/usr/bin/env bash
set -euo pipefail

# ./run_mafft.bash <input.fa> > <output.aln.fa>
#
# Thin wrapper around this farm's MAFFT install - /software/team301/mafft
# (the top-level wrapper script) is broken here (fails with bash "too
# many arguments" errors from unset internal path variables); the actual
# working binary is under mafft-7.525-with-extensions/scripts/mafft, and
# needs MAFFT_BINARIES pointed at its sibling binaries/ dir explicitly.
# Found and worked around 2026-09-22 while building the promoter-
# consensus step - see notes/plan.md.

MAFFT_BINARIES=/software/team301/mafft-7.525-with-extensions/binaries \
  /software/team301/mafft-7.525-with-extensions/scripts/mafft --auto --quiet "$1"
