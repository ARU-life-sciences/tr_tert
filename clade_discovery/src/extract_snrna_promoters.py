#!/usr/bin/env python3
"""Extracts the putative type-3 promoter region (upstream flank) of
each confident snRNA hit from a cmsearch run against Rfam's U1-U6/
U6atac/MRP/7SK models (see run_snrna_search_clade.bash). These are the
RNAPIII type-3-promoter-using snRNA families Fajkus used as anchors for
characterising a species/clade-specific USE (upstream sequence element)
consensus - this project's unknown TR is presumed to share the same
promoter architecture (that's the whole premise of the clade_discovery
effort: Hymenoptera/Lepidoptera TRs turned out to be this type, despite
sitting in Animalia), so a real USE consensus built from KNOWN
RNAPIII snRNAs in the same genome is the anchor for later finding the
TR itself by promoter architecture rather than (or in addition to) bare
template matching - see notes/plan.md's "Stage 0.5 candidate volume"
problem this was built to address.

Handles genomic strand correctly: for a '+' strand hit, "upstream" is
the lower-coordinate flank; for a '-' strand hit, the gene is read
right-to-left in genome coordinates, so "upstream" is the HIGHER-
coordinate flank, reverse-complemented back into transcript-relative
5'->3' orientation before writing.

Usage:
  extract_snrna_promoters.py <clade_name> <species> <tblout> [context_bp] [evalue_threshold]

Writes <outdir>/<species>_promoters.fa (one sequence per confident hit,
header: species|family|chrom|strand|evalue).
"""
import gzip
import sys
from pathlib import Path

CLADE_ROOT = Path(__file__).resolve().parent.parent
COMPLEMENT = str.maketrans("ACGT", "TGCA")
DEFAULT_CONTEXT = 200
DEFAULT_EVALUE = 1e-5


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def load_genome_path(clade_name, species):
    paths_file = CLADE_ROOT / "inputs" / f"{clade_name}_paths.txt"
    with open(paths_file) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == species:
                return parts[1]
    return None


def load_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs = {}
    name, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks).upper()
    return seqs


def parse_tblout(path, evalue_threshold):
    hits = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            chrom, family = parts[0], parts[2]
            seq_from, seq_to = int(parts[7]), int(parts[8])
            strand = parts[9]
            evalue = float(parts[15])
            if evalue < evalue_threshold:
                hits.append((chrom, family, seq_from, seq_to, strand, evalue))
    return hits


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    clade_name, species, tblout = sys.argv[1], sys.argv[2], sys.argv[3]
    context = int(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_CONTEXT
    evalue_threshold = float(sys.argv[5]) if len(sys.argv) > 5 else DEFAULT_EVALUE
    outdir = CLADE_ROOT / "outputs" / "snrna_search"
    outdir.mkdir(parents=True, exist_ok=True)

    genome_path = load_genome_path(clade_name, species)
    if genome_path is None:
        sys.exit(f"species {species!r} not found in {clade_name}_paths.txt")

    hits = parse_tblout(tblout, evalue_threshold)
    print(f"[extract_snrna_promoters] {species}: {len(hits)} confident hits "
          f"(e<{evalue_threshold})", file=sys.stderr)
    if not hits:
        return

    seqs = load_fasta(genome_path)
    out_fa = outdir / f"{species}_promoters.fa"
    n_written = 0
    with open(out_fa, "w") as fa:
        for chrom, family, seq_from, seq_to, strand, evalue in hits:
            if chrom not in seqs:
                continue
            chrom_seq = seqs[chrom]
            if strand == "+":
                gene_start = min(seq_from, seq_to)
                start = max(0, gene_start - context)
                promoter = chrom_seq[start:gene_start]
            else:
                gene_start = max(seq_from, seq_to)  # higher genomic coord = 5' end on - strand
                end = min(len(chrom_seq), gene_start + context)
                promoter = revcomp(chrom_seq[gene_start:end])
            if not promoter or "N" in promoter:
                continue
            header = f"{species}|{family}|{chrom}|{strand}|{evalue:.2e}"
            fa.write(f">{header}\n{promoter}\n")
            n_written += 1

    print(f"[extract_snrna_promoters] {species}: {n_written} promoter regions "
          f"-> {out_fa}", file=sys.stderr)


if __name__ == "__main__":
    main()
