#!/usr/bin/env bash
set -euo pipefail
# ----------------------------------------------------------------------
# run_sawfly_cm.bash
#
# The Tenthredinoidea TR covariance model, and a search of the sawflies
# still without a TR after model walking:
#   1. cmalign all 39 Tenthredinoidea TRs (10 published, 29 DToL) to the
#      seed CM (models/sawfly_tr.cm: 28 DToL TRs, MAFFT + RNAalifold),
#      cmbuild from that structural alignment, cmcalibrate
#   2. cmsearch (E <= 1e-3) of each target genome, hits with 300 bp
#      flanks -> outputs/sawfly_cm/<species>.hits.fa (as tr_hym), for
#      template_check.py
# ----------------------------------------------------------------------
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
I=/software/team301/infernal-1.1.4-linux-intel-gcc/binaries
H=/software/team301/hmmer-3.4
cd "$ROOT/models"
"$I/cmalign" --cpu 8 -o tenthredinoidea_all.sto sawfly_tr.cm ../../tr_hym_r5/models/Tenthredinoidea.fa > /dev/null
"$I/cmbuild" -F -n tenthredinoidea_tr tenthredinoidea_tr.cm tenthredinoidea_all.sto > /dev/null
"$I/cmcalibrate" --cpu 8 tenthredinoidea_tr.cm > tenthredinoidea_tr.cmcalibrate.log
"$I/cmpress" -F tenthredinoidea_tr.cm > /dev/null
out="$ROOT/outputs/sawfly_cm"; mkdir -p "$out"
SP="$ROOT/../tr_hym/species.tsv"
for sp in Aproceros_leucopoda Arge_ochropus Arge_ustulata Athalia_circularis Athalia_cordata Abia_aenea Urocerus_augur Allantus_rufocinctus Xiphydria_camelus; do
    f=$(awk -F'\t' -v s="$sp" '$1==s{print $2}' "$SP")
    t=$(mktemp -d "${TMPDIR:-/tmp}/scm.XXXXXX")
    gzip -dcf "$f" | tr -d '\000' > "$t/g.fa"
    "$H/easel/miniapps/esl-sfetch" --index "$t/g.fa" > /dev/null
    "$I/cmsearch" --cpu 8 -E 1e-3 --tblout "$t/cm.tbl" tenthredinoidea_tr.cm "$t/g.fa" > /dev/null
    awk '!/^#/ {print $1, "-", $3, "-", "-", "-", $8, $9, "-", "-", "-", "-", $16}' "$t/cm.tbl" > "$out/${sp}.nhmmer.tbl"
    { grep -v '^#' "$out/${sp}.nhmmer.tbl" || true; } | awk '{print $1, $7, $8, $3, $13}' | while read -r seq a b model e; do
        if (( a < b )); then s=$(( a > 300 ? a - 300 : 1 )); e2=$(( b + 300 )); else s=$(( a + 300 )); e2=$(( b > 300 ? b - 300 : 1 )); fi
        "$H/easel/miniapps/esl-sfetch" -n "${seq}:${a}-${b}|${model}|${e}" -c "${s}..${e2}" "$t/g.fa" "$seq" 2>/dev/null || true
    done > "$out/${sp}.hits.fa"
    rm -rf "$t"
    echo "[run_sawfly_cm] $sp: $(grep -c '>' "$out/${sp}.hits.fa") hits"
done
