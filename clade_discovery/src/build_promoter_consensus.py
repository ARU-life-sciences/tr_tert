#!/usr/bin/env python3
"""Computes a per-column consensus + conservation profile from a MAFFT
alignment of extracted snRNA promoter regions (extract_snrna_promoters.py
-> mafft), to characterise the type-3 promoter USE/TATA signature for a
species/clade - see notes/plan.md. Promoter sequences are extracted so
the transcription start site (gene start) sits at the 3' (right) end of
each sequence, so the proximal promoter elements (TATA box ~25-30bp
upstream, USE/PSE further upstream) should show up as conserved columns
concentrated toward the right end of the alignment.

Usage:
  build_promoter_consensus.py <alignment.fa> [min_frac=0.5]

Prints the per-column consensus (uppercase where conservation >=
min_frac, lowercase otherwise) and a compact conservation-fraction
track underneath, right-aligned so the 3' (TSS-proximal) end lines up.
"""
import sys
from collections import Counter


def load_fasta(path):
    seqs = []
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(chunks)))
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs.append((name, "".join(chunks)))
    return seqs


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    aln_path = sys.argv[1]
    min_frac = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5

    seqs = load_fasta(aln_path)
    n = len(seqs)
    aln_len = len(seqs[0][1])
    print(f"[build_promoter_consensus] {n} sequences, alignment length {aln_len}", file=sys.stderr)

    consensus = []
    fracs = []
    for col in range(aln_len):
        bases = [s[1][col] for s in seqs if s[1][col] != "-"]
        if not bases:
            consensus.append("-")
            fracs.append(0.0)
            continue
        counts = Counter(bases)
        base, count = counts.most_common(1)[0]
        frac = count / n  # fraction of ALL sequences (gaps count against), not just non-gap
        consensus.append(base if frac >= min_frac else base.lower())
        fracs.append(frac)

    cons_str = "".join(consensus)
    frac_track = "".join(str(int(f * 10)) if f < 1.0 else "*" for f in fracs)

    print(f"consensus (upper = >={min_frac*100:.0f}% of all {n} sequences agree):")
    print(cons_str)
    print("conservation (0-9 = tenths, * = 100%):")
    print(frac_track)

    # Flag the best-conserved stretch of >=6 columns with mean frac >= min_frac
    best = None
    for start in range(aln_len - 5):
        window = fracs[start:start + 8]
        mean_frac = sum(window) / len(window)
        if mean_frac >= min_frac and (best is None or mean_frac > best[0]):
            best = (mean_frac, start)
    if best:
        mean_frac, start = best
        print(f"\nbest-conserved 8bp window: position {start}-{start+8}, "
              f"mean conservation {mean_frac*100:.0f}%: {cons_str[start:start+8]}", file=sys.stderr)


if __name__ == "__main__":
    main()
