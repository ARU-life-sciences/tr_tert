#!/usr/bin/env python3
"""Synteny-anchored TR search, tested on pairs where both TRs are known.

TR sequence outside the template turns over too fast to align beyond close
relatives (Cerceris vs Mimumesa: no blastn similarity), but the TR's
genomic neighbourhood can stay collinear. So: map chunks of species A's
flanks around its TR to species B, take the anchors either side of A's TR
on B's best-matching sequence, and search only the interval between them
for B's template.

  1. chunks of CHUNK nt from A, up to FLANK nt either side of A's TR
     (skipping 1 kb next to it, which is the TR itself and its promoter)
  2. minimap2 -x asm20 of the chunks to B (primary hits, mapq >= 5)
  3. B's sequence with the most anchors; the anchors nearest A's TR on
     either side give B's interval (one side only: the nearest anchor
     +/- its distance from A's TR, doubled)
  4. template-like windows of B's unit in the interval (as covariation.py:
     len(unit) + 2 nt of its tandem, <= 3 copies), and whether one of them
     is B's known TR

Usage: synteny_test.py <genome_A> <A_TR seq:from-to> <genome_B> <unit_B> <B_TR seq:from-to> <outdir>
Prints one tab-separated line:
  anchors  B_seq  interval  interval_kb  B_TR_inside  windows  B_TR_window
"""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from covariation import read_fasta, rc, tandem_match, template_regex, MAX_COPIES  # noqa: E402

MINIMAP2 = str(Path(__file__).resolve().parents[2] / "bin" / "minimap2")
CHUNK = 2000
FLANK = 100_000
SKIP = 1000


def parse_locus(s):
    seq, rng = s.rsplit(":", 1)
    a, b = (int(x) for x in rng.split("-"))
    return seq, min(a, b), max(a, b)


def windows(seq, unit, offset):
    """Template-like windows of unit in seq (both strands), as genomic
    start positions (1-based) of B."""
    pat, w = template_regex(unit)
    out = []
    for strand, s in (("+", seq), ("-", rc(seq))):
        for m in pat.finditer(s):
            i = m.start()
            ctx = s[max(0, i - 3 * len(unit)): i + w + 3 * len(unit)]
            if tandem_match(ctx, unit) > MAX_COPIES * len(unit):
                continue
            pos = i if strand == "+" else len(s) - i - w
            out.append(offset + pos + 1)
    return sorted(set(out))


def main():
    gA, locA, gB, unitB, locB, outdir = sys.argv[1:7]
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    seqA, a0, a1 = parse_locus(locA)
    seqB_tr, b0, b1 = parse_locus(locB)
    A = dict(read_fasta(gA))[seqA]
    chunks = out / "chunks.fa"
    with open(chunks, "w") as fh:
        for start in range(max(0, a0 - FLANK), a0 - SKIP - CHUNK, CHUNK):
            fh.write(f">L{a0 - (start + CHUNK)}\n{A[start:start + CHUNK]}\n")
        for start in range(a1 + SKIP, min(len(A), a1 + FLANK) - CHUNK, CHUNK):
            fh.write(f">R{start - a1}\n{A[start:start + CHUNK]}\n")
    paf = subprocess.run([MINIMAP2, "-x", "asm20", "-t", "8", "--secondary=no", gB, chunks],
                         capture_output=True, text=True, check=True).stdout
    hits = {}
    for line in paf.splitlines():
        f = line.split("\t")
        if int(f[11]) < 5:
            continue
        name, tseq, ts, te = f[0], f[5], int(f[7]), int(f[8])
        if name not in hits or int(f[9]) > hits[name][3]:
            hits[name] = (tseq, ts, te, int(f[9]))
    if not hits:
        print("0\t-\t-\t-\tno_anchors\t-\t-")
        return
    by_seq = {}
    for name, (tseq, ts, te, _) in hits.items():
        by_seq.setdefault(tseq, []).append((name[0], int(name[1:]), ts, te))
    bseq, anchors = max(by_seq.items(), key=lambda x: len(x[1]))
    left = min((a for a in anchors if a[0] == "L"), key=lambda a: a[1], default=None)
    right = min((a for a in anchors if a[0] == "R"), key=lambda a: a[1], default=None)
    if left and right:
        lo, hi = min(left[2], right[2]), max(left[3], right[3])
    else:
        a = left or right
        lo, hi = a[2] - 2 * a[1] - CHUNK, a[3] + 2 * a[1] + CHUNK
    lo = max(0, lo)
    B = dict(read_fasta(gB))[bseq]
    hi = min(len(B), hi)
    inside = bseq == seqB_tr and lo <= b0 and b1 <= hi
    ws = windows(B[lo:hi], unitB, lo)
    in_tr = sum(1 for p in ws if b0 - 20 <= p <= b1 + 20)
    print(f"{len(anchors)}\t{bseq}\t{lo}-{hi}\t{(hi - lo) / 1000:.1f}\t"
          f"{'yes' if inside else 'no'}\t{len(ws)}\t{in_tr}")


if __name__ == "__main__":
    main()
