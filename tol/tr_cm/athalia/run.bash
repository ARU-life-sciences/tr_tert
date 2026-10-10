set -e
D=$(dirname "$0"); B=/software/team301/ncbi-blast-2.16.0+/bin; H=/software/team301/hmmer-3.4/src
$H/hmmbuild --dna "$D/ath.hmm" "$D/Athalia_rosae_tr.fa" > /dev/null
for sp in Athalia_circularis Athalia_cordata; do
  f=$(awk -F'\t' -v s=$sp '$1==s{print $2}' "$D/../../tr_hym/species.tsv")
  t=$(mktemp -d ${TMPDIR:-/tmp}/ath.XXXX); gzip -dcf "$f" | tr -d '\000' > $t/g.fa
  echo "== $sp blastn"; $B/blastn -query "$D/Athalia_rosae_tr.fa" -subject $t/g.fa -evalue 1e-5 -outfmt "6 sseqid sstart send pident length evalue" | head -5
  echo "== $sp nhmmer"; $H/nhmmer --dna -E 1e-3 --tblout $t/n.tbl "$D/ath.hmm" $t/g.fa > /dev/null; grep -v "^#" $t/n.tbl | awk '{print $1,$7,$8,$13}' | head -5
  rm -rf $t
done
