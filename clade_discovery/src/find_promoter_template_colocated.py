#!/usr/bin/env python3
"""Step 3, take 2: the actual fix for the candidate-volume problem that
motivated building the promoter-anchor machinery (see notes/plan.md).
Bare template-matching alone produced hundreds of thousands of
candidates per genome (short repeat units aren't specific enough on
their own - see predict_tr_candidate_loci.py's docstring). This inverts
the search: scan for the much rarer promoter-like signature first (a
near-match to the TATA-box-like consensus characterized from confirmed
snRNA genes via extract_snrna_promoters.py + build_promoter_consensus.py),
then only check the ~200bp window downstream of each promoter hit for a
template-region match - directly mirroring Fajkus's stated method:
"co-occurrence of both elements - the USE and the TR template region
(within 200 nt downstream of the USE) - enabled the identification of
TR candidates."

Usage:
  find_promoter_template_colocated.py <clade_name> <species> <repeat_unit> <promoter_consensus> [max_mismatches] [window_bp] [outdir]

  promoter_consensus: the short (~6-10bp) conserved core sequence found
  by build_promoter_consensus.py (e.g. "GTATATAT").
  max_mismatches: allowed Hamming-distance mismatches when scanning the
  genome for promoter-like near-matches (default 1).

Writes <outdir>/<species>_colocated_candidates.fa/.tsv - candidates
where a promoter-like match and a template-like match co-occur within
window_bp, in the correct relative orientation.
"""
import gzip
import re
import sys
from pathlib import Path

CLADE_ROOT = Path(__file__).resolve().parent.parent
COMPLEMENT = str.maketrans("ACGT", "TGCA")
DEFAULT_MAX_MM = 1
DEFAULT_WINDOW = 200
CONTEXT = 50  # extra flanking context around the final candidate pair


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def rotations(unit):
    return [unit[i:] + unit[:i] for i in range(len(unit))]


def all_variants(unit):
    variants = []
    seen = set()
    for v in rotations(unit):
        if v not in seen:
            variants.append(v)
            seen.add(v)
    for v in rotations(revcomp(unit)):
        if v not in seen:
            variants.append(v)
            seen.add(v)
    return variants


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


def build_fuzzy_regex(motif, max_mm):
    """Compiles a single alternation regex matching every sequence
    within `max_mm` Hamming distance of `motif` (e.g. for an 8bp motif
    and max_mm=1: 1 + 8*3 = 25 variants) - a pure-Python per-position
    sliding-window Hamming scan over a genome-scale sequence (hundreds
    of millions of bp) was ~40x too slow in practice (killed after 10+
    minutes CPU with no result on one genome); a compiled alternation
    regex runs the same search in C, comfortably fast enough."""
    bases = "ACGT"
    variants = {motif}
    if max_mm >= 1:
        for i in range(len(motif)):
            for b in bases:
                if b != motif[i]:
                    variants.add(motif[:i] + b + motif[i + 1:])
    if max_mm >= 2:
        for v in list(variants):
            for i in range(len(v)):
                for b in bases:
                    if b != v[i]:
                        variants.add(v[:i] + b + v[i + 1:])
    return re.compile("|".join(sorted(variants)))


def find_near_matches(seq, pattern):
    return [m.start() for m in pattern.finditer(seq) if "N" not in m.group()]


def has_template_match(window, variants):
    return any(v in window for v in variants)


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    clade_name, species, repeat_unit, promoter_consensus = sys.argv[1:5]
    repeat_unit = repeat_unit.upper()
    promoter_consensus = promoter_consensus.upper()
    max_mm = int(sys.argv[5]) if len(sys.argv) > 5 else DEFAULT_MAX_MM
    window = int(sys.argv[6]) if len(sys.argv) > 6 else DEFAULT_WINDOW
    outdir = Path(sys.argv[7]) if len(sys.argv) > 7 else CLADE_ROOT / "outputs" / "colocated_candidates"
    outdir.mkdir(parents=True, exist_ok=True)

    genome_path = load_genome_path(clade_name, species)
    if genome_path is None:
        sys.exit(f"species {species!r} not found in {clade_name}_paths.txt")

    print(f"[find_promoter_template_colocated] {species}: loading genome", file=sys.stderr)
    seqs = load_fasta(genome_path)
    variants = all_variants(repeat_unit)
    promoter_rc = revcomp(promoter_consensus)
    fwd_pattern = build_fuzzy_regex(promoter_consensus, max_mm)
    rev_pattern = build_fuzzy_regex(promoter_rc, max_mm)
    print(f"[find_promoter_template_colocated] {species}: fuzzy patterns built "
          f"({fwd_pattern.pattern.count('|') + 1} fwd variants)", file=sys.stderr)

    fa_out = outdir / f"{species}_colocated_candidates.fa"
    tsv_out = outdir / f"{species}_colocated_candidates.tsv"
    n_promoter_hits = 0
    n_candidates = 0

    with open(fa_out, "w") as fa, open(tsv_out, "w") as tsv:
        tsv.write("species\tchrom\tpromoter_pos\tstrand\ttemplate_pos\ttemplate_variant\n")
        for chrom, seq in seqs.items():
            # forward-strand promoter (transcript reads left-to-right,
            # template expected downstream = higher coordinate)
            for pos in find_near_matches(seq, fwd_pattern):
                n_promoter_hits += 1
                ds_start = pos + len(promoter_consensus)
                ds_end = min(len(seq), ds_start + window)
                downstream = seq[ds_start:ds_end]
                if not has_template_match(downstream, variants):
                    continue
                ctx_start = max(0, pos - CONTEXT)
                ctx_end = ds_end
                region = seq[ctx_start:ctx_end]
                header = f"{species}|{chrom}|{pos}|+"
                fa.write(f">{header}\n{region}\n")
                tsv.write(f"{species}\t{chrom}\t{pos}\t+\t{ds_start}\t{repeat_unit}\n")
                n_candidates += 1

            # reverse-strand promoter (transcript reads right-to-left in
            # genome coordinates, template expected at LOWER coordinate)
            for pos in find_near_matches(seq, rev_pattern):
                n_promoter_hits += 1
                ds_end = pos
                ds_start = max(0, ds_end - window)
                downstream = seq[ds_start:ds_end]
                if not has_template_match(downstream, variants):
                    continue
                ctx_start = ds_start
                ctx_end = min(len(seq), pos + len(promoter_rc) + CONTEXT)
                region = revcomp(seq[ctx_start:ctx_end])
                header = f"{species}|{chrom}|{pos}|-"
                fa.write(f">{header}\n{region}\n")
                tsv.write(f"{species}\t{chrom}\t{pos}\t-\t{ds_start}\t{repeat_unit}\n")
                n_candidates += 1

    print(f"[find_promoter_template_colocated] {species}: {n_promoter_hits} promoter-like "
          f"hits genome-wide, {n_candidates} with a co-located template match "
          f"-> {fa_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
