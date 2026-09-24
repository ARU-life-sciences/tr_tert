#!/usr/bin/env python3
"""Step 5 (structural validation) for the final, lineage-specific
candidate shortlist (see notes/plan.md's "Extended conservation check"
section). For each candidate, pulls the orthologous sequence from every
species it was confirmed conserved in (via the reciprocal blastn hit
coordinates already computed), aligns them with MAFFT, and runs
RNA.alifold (ViennaRNA's alignment-based consensus folding - installed
this session via `pip install --user ViennaRNA`, since no RNA structure
tool was available on this farm) to get a covariation-supported
consensus secondary structure, not just a single-sequence fold.

This is a real but partial structural check: standard nearest-neighbor
folding (what alifold does) cannot predict pseudoknots, so it can only
support/refute the presence of a plausible SIMPLE stem-loop (candidate
P1c-like element) - not the full PK+TBE+P1c architecture Fajkus
describes. Stated as a limitation, not glossed over.

Usage:
  structural_validation.py <clade_name> <final_candidates_file> <granarius_candidates_fa> <blast_tsv1> [blast_tsv2 ...]

  final_candidates_file: one query ID per line (species|chrom|pos|strand)
  blast_tsvN: reciprocal_filter_candidates.bash output files to pull
  orthologous coordinates from (subject species inferred from filename:
  <query_species>_vs_<subject_species>.tsv)

Writes <outdir>/<query_id_safe>.aln.fa and a summary TSV:
outputs/structural_validation/summary.tsv (candidate, n_species,
mean_covariation_pairs, has_plausible_stem_loop, mfe).
"""
import gzip
import re
import subprocess
import sys
from pathlib import Path

CLADE_ROOT = Path(__file__).resolve().parent.parent
COMPLEMENT = str.maketrans("ACGT", "TGCA")
MIN_STEM_PAIRS = 3  # a hairpin needs at least this many consecutive
                     # paired positions to count as a "plausible stem"


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


_genome_cache = {}


def load_fasta_cached(path):
    if path in _genome_cache:
        return _genome_cache[path]
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
    _genome_cache[path] = seqs
    return seqs


def load_granarius_candidates(fa_path):
    seqs = {}
    name, chunks = None, []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs


def parse_blast_hits(tsv_path):
    """Returns {query: (subject_chrom, sstart, send)} - best (first,
    already filtered upstream to single-hit) row per query."""
    hits = {}
    with open(tsv_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            q, sseqid, sstart, send = parts[0], parts[1], int(parts[8]), int(parts[9])
            hits[q] = (sseqid, sstart, send)
    return hits


def extract_ortholog(genome_path, chrom, sstart, send, pad=20):
    seqs = load_fasta_cached(genome_path)
    if chrom not in seqs:
        return None
    seq = seqs[chrom]
    lo, hi = min(sstart, send), max(sstart, send)
    lo = max(0, lo - 1 - pad)
    hi = min(len(seq), hi + pad)
    region = seq[lo:hi]
    if sstart > send:  # hit was on the minus strand relative to query
        region = revcomp(region)
    return region


def run_mafft(fa_path, aln_path):
    result = subprocess.run(
        ["bash", str(CLADE_ROOT / "src" / "run_mafft.bash"), str(fa_path)],
        capture_output=True, text=True,
    )
    aln_path.write_text(result.stdout)
    return result.returncode == 0 and result.stdout.strip() != ""


def analyse_structure(aln_path):
    import RNA
    names, seqs = [], []
    name, chunks = None, []
    with open(aln_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs.append("".join(chunks))
                name = line[1:].strip()
                names.append(name)
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs.append("".join(chunks))
    if len(seqs) < 2:
        return None
    structure, mfe, *_ = RNA.alifold(seqs)
    # count consecutive-paired-position runs as candidate stems
    pairs = []
    stack = []
    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            if stack:
                pairs.append((stack.pop(), i))
    longest_run = 0
    if pairs:
        pairs.sort()
        run = 1
        for i in range(1, len(pairs)):
            if pairs[i][0] == pairs[i - 1][0] + 1 and pairs[i][1] == pairs[i - 1][1] - 1:
                run += 1
            else:
                longest_run = max(longest_run, run)
                run = 1
        longest_run = max(longest_run, run)
    return dict(n_seqs=len(seqs), structure=structure, mfe=mfe,
                n_pairs=len(pairs), longest_stem_run=longest_run)


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    clade_name = sys.argv[1]
    candidates_file = sys.argv[2]
    granarius_fa = sys.argv[3]
    blast_tsvs = sys.argv[4:]

    candidates = [l.strip() for l in open(candidates_file) if l.strip()]
    granarius_seqs = load_granarius_candidates(granarius_fa)

    # subject_species -> {query: (chrom, sstart, send)}
    per_species_hits = {}
    for tsv in blast_tsvs:
        m = re.match(r".*_vs_(.+)\.tsv$", Path(tsv).name)
        if not m:
            print(f"Skipping {tsv}: can't parse subject species from filename", file=sys.stderr)
            continue
        subject_sp = m.group(1)
        per_species_hits[subject_sp] = parse_blast_hits(tsv)

    outdir = CLADE_ROOT / "outputs" / "structural_validation"
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "summary.tsv"

    with open(summary_path, "w") as summary:
        summary.write("candidate\tn_species\tmfe\tn_covarying_pairs\tlongest_stem_run\tplausible_stem_loop\n")
        for i, cand in enumerate(candidates, 1):
            safe_name = cand.replace("|", "_")
            print(f"[structural_validation] ({i}/{len(candidates)}) {cand}", file=sys.stderr)
            fa_path = outdir / f"{safe_name}.fa"
            with open(fa_path, "w") as fa:
                fa.write(f">Aphodius_granarius\n{granarius_seqs[cand]}\n")
                for subject_sp, hits in per_species_hits.items():
                    if cand not in hits:
                        continue
                    chrom, sstart, send = hits[cand]
                    genome_path = load_genome_path(clade_name, subject_sp)
                    ortholog = extract_ortholog(genome_path, chrom, sstart, send)
                    if ortholog:
                        fa.write(f">{subject_sp}\n{ortholog}\n")

            n_seqs = sum(1 for _ in open(fa_path) if _.startswith(">"))
            if n_seqs < 2:
                summary.write(f"{cand}\t{n_seqs}\tNA\tNA\tNA\tFALSE\n")
                continue

            aln_path = outdir / f"{safe_name}.aln.fa"
            ok = run_mafft(fa_path, aln_path)
            if not ok:
                summary.write(f"{cand}\t{n_seqs}\tNA\tNA\tNA\tFALSE\n")
                continue

            result = analyse_structure(aln_path)
            if result is None:
                summary.write(f"{cand}\t{n_seqs}\tNA\tNA\tNA\tFALSE\n")
                continue
            plausible = result["longest_stem_run"] >= MIN_STEM_PAIRS
            summary.write(f"{cand}\t{result['n_seqs']}\t{result['mfe']:.2f}\t"
                           f"{result['n_pairs']}\t{result['longest_stem_run']}\t{plausible}\n")

    print(f"[structural_validation] done -> {summary_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
