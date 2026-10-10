#!/usr/bin/env bash
# NCBI Datasets: genome FASTA for each assembly in inputs/assemblies.tsv
# (plus protein FASTA and GFF3 where annotated) -> genomes/<species>/
set -uo pipefail
cd "$(dirname "$0")/.."
tail -n +2 inputs/assemblies.tsv | while IFS=$'\t' read -r sp acc level note; do
    d="genomes/$sp"; mkdir -p "$d"
    [[ -s "$d/$acc.fa.gz" ]] && continue
    types="GENOME_FASTA"; [[ $acc == GCF_* ]] && types="GENOME_FASTA,PROT_FASTA,GENOME_GFF"
    curl -sL -m 7200 -o "$d/$acc.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=$types"
    unzip -q -o "$d/$acc.zip" -d "$d/tmp" && \
      find "$d/tmp" -name "*_genomic.fna" -exec cat {} + | gzip > "$d/$acc.fa.gz" && \
      { find "$d/tmp" -name "protein.faa" -exec cp {} "$d/protein.faa" \; ; find "$d/tmp" -name "genomic.gff" -exec cp {} "$d/genomic.gff" \; ; true; }
    rm -rf "$d/tmp" "$d/$acc.zip"
    echo "[download] $sp $acc $(du -h "$d/$acc.fa.gz" | cut -f1)"
done
