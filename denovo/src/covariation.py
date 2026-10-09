#!/usr/bin/env python3
"""De novo TR candidates by template co-variation across a repeat shift.

A TR's template is ~1.2-2 copies of its species' telomeric repeat (C-rich
strand, TR sense). Where the repeat changed between two related species,
the orthologous TR locus must have a template matching each species' own
repeat: the shift unit in one, the other unit in the other. Template-like
sequence that is there by chance has no reason to change with the
telomere, so requiring the change filters it out, whatever the TR's
biogenesis (snRNA-type, H/ACA) or promoter.

Steps:
  1. candidates in species A: windows of len(unit_A) + MIN_EXTRA nt
     exactly matching the tandem of unit_A (any phase, either strand), not
     inside a longer array (the telomere and interstitial telomeric
     sequence: <= MAX_COPIES copies), with FLANK nt either side
  2. blastn of the candidates against species B's genome (E <= 1e-10)
  3. for each hit spanning the template window, the aligned window in B:
     a co-varying candidate if B's window matches unit_B (>= len(unit_B) +
     MIN_EXTRA) and not unit_A (< len(unit_A) + MIN_EXTRA)

If unit_A == unit_B (no shift between the two species), step 3 instead
keeps candidates whose aligned window in B still matches the unit: a TR
keeps its template, so loci conserved with the template intact across
several species of a group are TR candidates (conservation mode; weaker
than co-variation, as any conserved template-like sequence passes).

Usage (on a node with ~16 GB; makeblastdb on B):
  covariation.py <genome_A> <unit_A> <genome_B> <unit_B> <outdir> [known_tr.tsv]

known_tr.tsv: optional seq<TAB>start<TAB>end of A's known TR, to report
whether (and at what rank) the known locus is recovered.

Outputs in outdir: covarying.tsv (ranked by blast bitscore) and
summary.txt; the candidates, blast hits and B's database are deleted at the
end (tens of GB for short units).
"""
import gzip
import re
import subprocess
import sys
from pathlib import Path

BLAST = "/software/team301/ncbi-blast-2.16.0+/bin"
MIN_EXTRA = 2
MAX_COPIES = 3
FLANK = 150
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    op = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name:
        yield name, "".join(chunks).upper()


def tandem_match(seq, unit):
    n, best = len(unit), 0
    for start in range(len(seq)):
        for phase in range(n):
            if seq[start] != unit[phase]:
                continue
            k = 0
            while start + k < len(seq) and seq[start + k] == unit[(phase + k) % n]:
                k += 1
            best = max(best, k)
    return best


def template_regex(unit):
    """All len(unit)+MIN_EXTRA windows of the unit's tandem."""
    w = len(unit) + MIN_EXTRA
    t = unit * (w // len(unit) + 2)
    words = {t[p:p + w] for p in range(len(unit))}
    return re.compile("(?=(" + "|".join(sorted(words)) + "))"), w


def candidates(genome, unit, out_fa):
    """Template-like windows of unit in genome, both strands, outside
    arrays; writes them with flanks (TR sense: the C-rich strand)."""
    pat, w = template_regex(unit)
    n = 0
    with open(out_fa, "w") as fh:
        for name, seq in read_fasta(genome):
            for strand, s in (("+", seq), ("-", rc(seq))):
                last = -10**9
                for m in pat.finditer(s):
                    i = m.start()
                    if i - last < w:
                        last = i
                        continue
                    last = i
                    ctx = s[max(0, i - 3 * len(unit)): i + w + 3 * len(unit)]
                    if tandem_match(ctx, unit) > MAX_COPIES * len(unit):
                        continue
                    a, b = max(0, i - FLANK), i + w + FLANK
                    L = len(s)
                    # genomic coordinates (1-based, + strand) of the window
                    g0, g1 = (i + 1, i + w) if strand == "+" else (L - i - w + 1, L - i)
                    fh.write(f">{name}:{g0}-{g1}:{strand}|{i - a}\n{s[a:b]}\n")
                    n += 1
    return n


def main():
    gA, uA, gB, uB, outdir = sys.argv[1:6]
    known = sys.argv[6] if len(sys.argv) > 6 else None
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    n = candidates(gA, uA, out / "candidates.fa")
    # B as a blast database
    dbfa = out / "B.fa"
    with open(dbfa, "w") as fh:
        for name, seq in read_fasta(gB):
            fh.write(f">{name}\n{seq}\n")
    subprocess.run([f"{BLAST}/makeblastdb", "-in", dbfa, "-dbtype", "nucl", "-out", out / "B"],
                   check=True, stdout=subprocess.DEVNULL)
    subprocess.run([f"{BLAST}/blastn", "-query", out / "candidates.fa", "-db", out / "B",
                    "-evalue", "1e-10", "-max_target_seqs", "5", "-num_threads", "8",
                    "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qseq sseq",
                    "-out", out / "blast.tsv"], check=True)
    wA, wB = len(uA) + MIN_EXTRA, len(uB) + MIN_EXTRA
    cov = {}
    for line in open(out / "blast.tsv"):
        f = line.rstrip("\n").split("\t")
        q, off = f[0].rsplit("|", 1)
        off = int(off)
        qs, qe = int(f[4]), int(f[5])
        if not (qs <= off + 1 and qe >= off + wA):
            continue
        # walk the alignment to the window, then take B's aligned bases
        # plus a margin either side for a template of a different length
        qseq, sseq = f[10], f[11]
        qpos, start_col, end_col = qs - 1, None, None
        for col, c in enumerate(qseq):
            if c != "-":
                if qpos == off and start_col is None:
                    start_col = col
                if qpos == off + wA - 1:
                    end_col = col
                qpos += 1
        if start_col is None or end_col is None:
            continue
        lo, hi = max(0, start_col - 8), min(len(sseq), end_col + 9)
        bwin = sseq[lo:hi].replace("-", "")
        mB, mA = tandem_match(bwin, uB), tandem_match(bwin, uA)
        # same unit in both species: conservation mode, the template must be
        # retained; different units: co-variation, it must have changed
        keep = mB >= wB if uA == uB else (mB >= wB and mA < wA)
        if keep:
            key = q
            score = float(f[9])
            if key not in cov or score > cov[key][0]:
                cov[key] = (score, f[1], f[6], f[7], f[2], bwin, mB)
    ranked = sorted(cov.items(), key=lambda x: -x[1][0])
    with open(out / "covarying.tsv", "w") as fh:
        print("rank\tcandidate_A\tbitscore\tseq_B\tstart_B\tend_B\tpident\twindow_B\tunit_B_match", file=fh)
        for r, (q, v) in enumerate(ranked, 1):
            print(r, q, *v, sep="\t", file=fh)
    blasted = len({l.split("\t")[0] for l in open(out / "blast.tsv")})
    lines = [f"A unit {uA}: {n} template-like candidates; {blasted} with a blastn hit in B",
             f"B unit {uB}: {len(ranked)} co-varying candidates"]
    if known:
        kseq, ks, ke = open(known).read().split()[:3]
        ks, ke = sorted((int(ks), int(ke)))
        hit = [r for r, (q, _) in enumerate(ranked, 1)
               if q.split(":")[0] == kseq and ks - 300 <= int(q.split(":")[1].split("-")[0]) <= ke + 300]
        lines.append(f"known TR {kseq}:{ks}-{ke}: " + (f"recovered at rank {hit[0]}" if hit else "not recovered"))
    (out / "summary.txt").write_text("\n".join(lines) + "\n")
    # B's blast database and the raw candidates/hits are regenerable and
    # large (tens of GB for short units): keep only the results
    for f in out.glob("B.*"):
        f.unlink()
    for name in ("candidates.fa", "blast.tsv"):
        (out / name).unlink(missing_ok=True)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
