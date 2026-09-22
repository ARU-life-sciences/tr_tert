#!/usr/bin/env python3
"""Step 3 of the clade_discovery method (see notes/plan.md): predicts
candidate TR-like loci in a genome from a confirmed telomere repeat
unit, mirroring Fajkus's approach (Allium/Bombus): a genuine TR
template region is a near-single-copy genomic match to the repeat
(any rotation, either strand) - NOT part of the long tandem array that
IS the telomere/interstitial-repeat itself.

Specificity note (why this isn't just "search for the repeat unit",
and why a FIXED "exactly 2 copies" rule - the first version of this
script used exactly that, copied directly from the plant/Bombus
precedent - turned out to be wrong): those precedents used 11-13bp
repeat units, where even 2 tandem copies (22-26bp of exact match) is
genuinely rare by chance. This project's units have been much shorter
(this clade's confirmed motif is 5bp), and specificity scales as
4^(matched length) - a fixed "2 copies" rule for a 5bp unit only
requires 10bp of exact match, nowhere near enough: running this on a
real ~300Mb beetle genome (Aphodius_granarius) produced over 1 MILLION
"candidates" - almost entirely short-tandem/microsatellite noise, not
usable input for the reciprocal-blastn filtering step. Real genomes
also aren't uniformly random at this scale (satellite/low-complexity
DNA inflates short-tandem-repeat frequency well above a naive chance
estimate), which independently pushes hit counts up further.

Fix: scale the required copy count to the unit length so the total
matched span stays in a genuinely rare/specific range (~15bp,
MIN_SPECIFIC_SPAN below) regardless of how short the underlying unit
is, while still respecting the plant/Bombus principle of excluding the
real telomeric/interstitial repeat array itself - allow a narrow window
of MIN_COPIES to MIN_COPIES+1 tandem copies, not an open-ended "2 or
more". For a 5bp unit this means requiring 3-4 tandem copies (15-20bp),
not 2.

Usage:
  predict_tr_candidate_loci.py <clade_name> <species> <repeat_unit> [context_bp] [outdir]

Writes <outdir>/<species>_candidates.fa (context-extended candidate
sequences, one per locus) and <outdir>/<species>_candidates.tsv
(chrom, position, variant, orientation).
"""
import gzip
import re
import sys
from pathlib import Path

CLADE_ROOT = Path(__file__).resolve().parent.parent
COMPLEMENT = str.maketrans("ACGT", "TGCA")
DEFAULT_CONTEXT = 300
MERGE_DISTANCE = 20  # candidate hits from different rotations within this
                      # many bp of each other are treated as the same
                      # underlying locus (detected via an overlapping
                      # reading frame), not two independent candidates


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def rotations(unit):
    return [unit[i:] + unit[:i] for i in range(len(unit))]


def all_variants(unit):
    variants = []
    seen = set()
    for v in rotations(unit):
        if v not in seen:
            variants.append((v, "fwd"))
            seen.add(v)
    for v in rotations(revcomp(unit)):
        if v not in seen:
            variants.append((v, "rev"))
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


MIN_SPECIFIC_SPAN = 15  # bp - see module docstring for why this scales
                         # with unit length rather than being a fixed
                         # copy count


def copy_range(unit):
    min_copies = max(2, -(-MIN_SPECIFIC_SPAN // len(unit)))  # ceil div
    return min_copies, min_copies + 1


def find_hits(seq, unit):
    """Tandem duplications of `unit` whose total copy count falls in a
    narrow window (copy_range) - long enough to be specific for short
    units, short enough to still exclude the real telomeric/interstitial
    repeat array itself."""
    lo, hi = copy_range(unit)
    pattern = re.compile(f"(?<!{unit}){unit}{{{lo},{hi}}}(?!{unit})")
    return [m.start() for m in pattern.finditer(seq)]


def merge_hits(hits):
    """hits: list of (pos, variant, orientation). Merge hits within
    MERGE_DISTANCE bp into one candidate (keep the first seen)."""
    hits = sorted(hits, key=lambda h: h[0])
    merged = []
    for pos, variant, orient in hits:
        if merged and pos - merged[-1][0] <= MERGE_DISTANCE:
            continue
        merged.append((pos, variant, orient))
    return merged


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    clade_name, species, repeat_unit = sys.argv[1], sys.argv[2], sys.argv[3].upper()
    context = int(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_CONTEXT
    outdir = Path(sys.argv[5]) if len(sys.argv) > 5 else CLADE_ROOT / "outputs" / "candidate_loci"
    outdir.mkdir(parents=True, exist_ok=True)

    genome_path = load_genome_path(clade_name, species)
    if genome_path is None:
        sys.exit(f"species {species!r} not found in {clade_name}_paths.txt")

    print(f"[predict_tr_candidate_loci] {species}: loading genome", file=sys.stderr)
    seqs = load_fasta(genome_path)
    variants = all_variants(repeat_unit)
    print(f"[predict_tr_candidate_loci] {species}: {len(variants)} rotation/strand "
          f"variants of {repeat_unit}", file=sys.stderr)

    fa_out = outdir / f"{species}_candidates.fa"
    tsv_out = outdir / f"{species}_candidates.tsv"
    n_total = 0
    with open(fa_out, "w") as fa, open(tsv_out, "w") as tsv:
        tsv.write("species\tchrom\tposition\tvariant\torientation\n")
        for chrom, seq in seqs.items():
            all_hits = []
            for variant, orient in variants:
                for pos in find_hits(seq, variant):
                    all_hits.append((pos, variant, orient))
            merged = merge_hits(all_hits)
            for pos, variant, orient in merged:
                start = max(0, pos - context)
                _, hi_copies = copy_range(variant)
                end = min(len(seq), pos + len(variant) * hi_copies + context)
                window = seq[start:end]
                if "N" in window:
                    continue  # skip candidates sitting in/near an assembly gap
                header = f"{species}|{chrom}|{pos}|{variant}|{orient}"
                fa.write(f">{header}\n{window}\n")
                tsv.write(f"{species}\t{chrom}\t{pos}\t{variant}\t{orient}\n")
                n_total += 1

    print(f"[predict_tr_candidate_loci] {species}: {n_total} candidate loci "
          f"-> {fa_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
