#!/usr/bin/env python3
"""Are the extra TR-like copies retrocopies?

For each TR-like locus in a genome (hits.fa from a clade search, E <=
1e-5, one per 2 kb), with 1 kb flanks in the TR's orientation:
  polyA      a run of >= 10 A (1 mismatch allowed) starting within 60 nt
             after the copy's 3' end: the reverse-transcribed RNA's tail
  tsd        the longest exact direct repeat (>= 8 nt) with one copy within
             60 nt before the copy and one within 80 nt after it (after a
             poly-A, if any): a target-site duplication
  template   the copy carries a template for the species' tidk unit
  promoter   a TATA-like box (TATAWA/TATAAA, 1 mismatch) 20-40 nt before
             the copy's 5' end (the type-3 Pol III promoter element; it is
             upstream of the transcript, so not copied by retroposition)
Background: the same features at 20 random positions per genome (same
sequences, random orientation).

Usage: copies.py <species list: species fasta unit hits.fa> <out.tsv>
Output rows: species kind(copy/random) locus evalue template polyA tsd_len promoter
"""
import gzip
import random
import re
import sys

FLANK = 1000
COMP = str.maketrans("ACGTN", "TGCAN")
TATA = re.compile(r"(?=(TATA[AT]A|TATAAA|[ACG]ATA[AT]A|T[ACG]TA[AT]A|TA[ACG]A[AT]A|TAT[ACG][AT]A))")


def rc(s):
    return s.translate(COMP)[::-1]


def read_genome(path):
    op = gzip.open if path.endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip().replace("\x00", ""))
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def tandem(seq, unit):
    n, best = len(unit), 0
    for i in range(len(seq)):
        for p in range(n):
            if seq[i] != unit[p]:
                continue
            k = 0
            while i + k < len(seq) and seq[i + k] == unit[(p + k) % n]:
                k += 1
            best = max(best, k)
    return best


def polya_end(tail):
    """End index of a >= 10 nt A run (1 mismatch) starting within 60 nt, or None."""
    for start in range(0, 60):
        mism, k = 0, 0
        while start + k < len(tail):
            if tail[start + k] != "A":
                mism += 1
                if mism > 1:
                    break
            k += 1
        if k - mism >= 10:
            return start + k
    return None


def tsd(left, right, kmin=8):
    best = 0
    for k in range(kmin, 25):
        found = any(left[i:i + k] in right for i in range(len(left) - k + 1))
        if not found:
            break
        best = k
    return best


def features(seq, a, b, strand, unit):
    """seq: chromosome; a..b (0-based, b exclusive) the copy on + strand."""
    lo, hi = max(0, a - FLANK), min(len(seq), b + FLANK)
    region = seq[lo:hi]
    s, e = a - lo, b - lo
    if strand == "-":
        region = rc(region)
        s, e = len(region) - e, len(region) - s
    core = region[s:e]
    tail = region[e:e + 200]
    pa = polya_end(tail)
    right = region[e + (pa or 0):e + (pa or 0) + 80]
    left = region[max(0, s - 60):s]
    prom = bool(TATA.search(region[max(0, s - 40):max(0, s - 20)]))
    tmpl = bool(unit) and tandem(region[max(0, s - 300):e + 300], unit) >= len(unit) + 2
    return tmpl, pa is not None, tsd(left, right), prom


def main():
    species_list, out = sys.argv[1:3]
    rows = []
    random.seed(1)
    for line in open(species_list):
        if line.startswith("species"):
            continue
        sp, fasta, unit, hits = line.rstrip("\n").split("\t")
        loci = {}
        for h in (l[1:].strip() for l in open(hits) if l.startswith(">")):
            loc, model, e = h.split("|")
            seqname, rng = loc.rsplit(":", 1)
            x, y = (int(v) for v in rng.split("-"))
            if float(e) > 1e-5:
                continue
            key = (seqname, min(x, y) // 2000)
            if key not in loci or float(e) < loci[key][3]:
                loci[key] = (seqname, min(x, y) - 1, max(x, y), float(e), "+" if x < y else "-")
        if not loci:
            continue
        G = read_genome(fasta)
        for seqname, a, b, e, strand in loci.values():
            t, p, d, pr = features(G[seqname], a, b, strand, unit)
            rows.append((sp, "copy", f"{seqname}:{a}-{b}{strand}", f"{e:.1e}", t, p, d, pr))
        names = [n for n in G if len(G[n]) > 1_000_000]
        for _ in range(20):
            n = random.choice(names)
            a = random.randrange(FLANK, len(G[n]) - FLANK - 200)
            t, p, d, pr = features(G[n], a, a + 150, random.choice("+-"), unit)
            rows.append((sp, "random", f"{n}:{a}", "-", t, p, d, pr))
        del G
        print(f"[copies] {sp}: {len(loci)} copies", file=sys.stderr)
    with open(out, "w") as fh:
        print("species\tkind\tlocus\tevalue\ttemplate\tpolyA\ttsd_len\tpromoter", file=fh)
        for r in rows:
            print("\t".join(str(x) for x in r), file=fh)


if __name__ == "__main__":
    main()
