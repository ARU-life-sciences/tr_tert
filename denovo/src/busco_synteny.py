#!/usr/bin/env python3
"""Synteny-anchored TR search with BUSCO genes as anchors.

Nucleotide anchors (synteny_test.py) run out between families: intergenic
DNA diverges. Single-copy orthologous genes keep their identity, and ToL's
BUSCO pipeline has already placed ~5,800 Hymenoptera BUSCOs in every
released assembly. So:
  1. the K complete BUSCOs nearest species A's TR on either side
  2. the same BUSCOs in species B; the B sequence carrying most of them
  3. B's interval: between the nearest A-left and A-right BUSCOs found on
     that sequence (one side only: that gene +/- twice its distance in A)
  4. in the interval: nhmmer of a TR profile HMM and, optionally, cmsearch
     of a TR covariance model (E-values scale with the interval, not the
     genome); every hit's template is checked against B's unit

BUSCO tables are under <species dir>/assembly/release/reference/analysis/
pipelines/busco/*/results/busco/<lineage>/*.full_table.tsv.gz, in INSDC
sequence names, mapped to the curated SUPER_* names with the release's
name_mapping.tsv.

Usage: busco_synteny.py <fasta_A> <A_TR seq:from-to> <fasta_B> <unit_B> <B_TR seq:from-to | none>
                        <outdir> <tr.hmm> [tr.cm] [lineage=hymenoptera_odb12]
Prints one tab-separated line:
  anchors  B_seq  interval  interval_kb  B_TR_inside  best_hit  best_E  best_template  n_hits_with_template
"""
import glob
import gzip
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from covariation import read_fasta, rc, tandem_match  # noqa: E402

NHMMER = "/software/team301/hmmer-3.4/src/nhmmer"
CMSEARCH = "/software/team301/infernal-1.1.4-linux-intel-gcc/binaries/cmsearch"
K = 3
MAX_GAP = 1_000_000  # anchors further apart than this on B are separate blocks
HERE = Path(__file__).resolve().parent.parent
MINIPROT = str(HERE.parent / "bin" / "miniprot")
REF_PROTEINS = str(HERE / "inputs" / "busco_ref" / "hymenoptera_odb12_xiphydria.faa")
CACHE = HERE / "outputs" / "busco_tables"


def species_dir(fasta):
    return Path(fasta.split("/assembly/")[0])


def miniprot_table(fasta):
    """BUSCO id -> (seq, start, end) from miniprot of the reference BUSCO
    proteins (REF_PROTEINS: Xiphydria's hymenoptera_odb12 single-copy set)
    against an assembly with no BUSCO table of its own (e.g. a recent
    release). A protein's position is its best alignment, kept only if the
    second best scores < 90% of it (single copy). Cached per assembly."""
    tolid = fasta.split("/curated/")[1].split("/")[0]
    cache = CACHE / f"{tolid}.tsv"
    if not cache.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        paf = subprocess.run([MINIPROT, "-t", "8", "-I", "--outn=2", fasta, REF_PROTEINS],
                             capture_output=True, text=True, check=True).stdout
        best = {}
        for line in paf.splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 12:
                continue
            score = int(next(t[5:] for t in f[12:] if t.startswith("AS:i:")))
            best.setdefault(f[0], []).append((score, f[5], int(f[7]), int(f[8]),
                                             (int(f[3]) - int(f[2])) / int(f[1])))
        with open(cache, "w") as fh:
            for bid, hits in best.items():
                hits.sort(reverse=True)
                if hits[0][4] < 0.5 or (len(hits) > 1 and hits[1][0] >= 0.9 * hits[0][0]):
                    continue
                print(bid, *hits[0][1:4], sep="\t", file=fh)
    out = {}
    for line in open(cache):
        bid, seq, a, b = line.rstrip("\n").split("\t")
        out[bid] = (seq, int(a), int(b))
    return out


def busco_table(fasta, lineage):
    """BUSCO id -> (curated seq, start, end) for complete single-copy BUSCOs;
    from ToL's BUSCO run, or miniprot_table() if there isn't one."""
    ref = species_dir(fasta) / "assembly" / "release" / "reference"
    tables = glob.glob(str(ref / "analysis" / "pipelines" / "busco" / "*" / "results" / "busco"
                           / lineage / "*.full_table.tsv.gz"))
    maps = glob.glob(str(ref / "insdc" / "*.name_mapping.tsv"))
    if not tables or not maps:
        return miniprot_table(fasta) if lineage == "hymenoptera_odb12" else None
    names = {}
    for line in open(maps[0]):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            names[f[0]] = f[1]
    out = {}
    with gzip.open(sorted(tables)[-1], "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[1] != "Complete":
                continue
            a, b = int(f[3]), int(f[4])
            out[f[0]] = (names.get(f[2], f[2]), min(a, b), max(a, b))
    return out


def parse_locus(s):
    seq, rng = s.rsplit(":", 1)
    a, b = (int(x) for x in rng.split("-"))
    return seq, min(a, b), max(a, b)


def search(B_seq, lo, hi, seqname, unit, hmm, cm, out):
    """(best hit, its E, its template ok, number of hits with a template);
    hits in genomic coordinates on seqname."""
    iv = out / "interval.fa"
    iv.write_text(f">{seqname}_{lo}\n{B_seq[lo:hi]}\n")
    hits = []
    tbl = out / "interval.nhmmer.tbl"
    subprocess.run([NHMMER, "--dna", "-E", "10", "--tblout", tbl, hmm, iv],
                   check=True, stdout=subprocess.DEVNULL)
    for l in open(tbl):
        if not l.startswith("#"):
            f = l.split()
            hits.append(("hmm", int(f[6]), int(f[7]), float(f[12])))
    if cm:
        tbl = out / "interval.cm.tbl"
        subprocess.run([CMSEARCH, "-E", "10", "--tblout", tbl, cm, iv],
                       check=True, stdout=subprocess.DEVNULL)
        for l in open(tbl):
            if not l.startswith("#"):
                f = l.split()
                hits.append(("cm", int(f[7]), int(f[8]), float(f[15])))
    scored = []
    for kind, a, b, e in hits:
        g0, g1 = lo + min(a, b), lo + max(a, b)
        region = B_seq[max(0, g0 - 300):g1 + 300]
        if a > b:
            region = rc(region)
        ok = tandem_match(region, unit) >= len(unit) + 2
        scored.append((e, kind, g0, g1, ok))
    scored.sort()
    with open(out / "interval_hits.tsv", "w") as fh:
        for e, kind, g0, g1, ok in scored:
            print(kind, f"{seqname}:{g0}-{g1}", f"{e:.2e}", "template" if ok else "-", sep="\t", file=fh)
    if not scored:
        return "-", "-", "-", 0
    e, kind, g0, g1, ok = scored[0]
    return (f"{seqname}:{g0}-{g1}({kind})", f"{e:.1e}", "yes" if ok else "no",
            sum(1 for x in scored if x[4]))


def main():
    gA, locA, gB, unitB, locB, outdir, hmm = sys.argv[1:8]
    cm = sys.argv[8] if len(sys.argv) > 8 and sys.argv[8] != "-" else None
    lineage = sys.argv[9] if len(sys.argv) > 9 else "hymenoptera_odb12"
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    tA, tB = busco_table(gA, lineage), busco_table(gB, lineage)
    if not tA or not tB:
        print("0\t-\t-\t-\tno_busco_table\t-\t-\t-\t0")
        return
    seqA, a0, a1 = parse_locus(locA)
    onA = [(bid, s, e) for bid, (q, s, e) in tA.items() if q == seqA]
    left = sorted([x for x in onA if x[2] < a0], key=lambda x: a0 - x[2])[:K]
    right = sorted([x for x in onA if x[1] > a1], key=lambda x: x[1] - a1)[:K]
    anchors = [(side, bid, tB[bid]) for side, xs in (("L", left), ("R", right))
               for bid, _, _ in xs if bid in tB]
    if not anchors:
        print("0\t-\t-\t-\tno_anchors\t-\t-\t-\t0")
        return
    # the densest collinear block: anchors on one B sequence, split where
    # consecutive anchors are > MAX_GAP apart (rearrangements otherwise
    # stretch the interval to megabases)
    blocks = []
    for q in {a[2][0] for a in anchors}:
        on = sorted((a for a in anchors if a[2][0] == q), key=lambda a: a[2][1])
        cur = [on[0]]
        for a in on[1:]:
            if a[2][1] - cur[-1][2][2] > MAX_GAP:
                blocks.append(cur)
                cur = []
            cur.append(a)
        blocks.append(cur)
    on = max(blocks, key=len)
    bseq = on[0][2][0]
    in_block = {a[1] for a in on}
    # nearest-in-A anchor on each side within the block
    near = {}
    for side, xs in (("L", left), ("R", right)):
        for bid, s0, e0 in xs:
            if bid in in_block:
                near[side] = (tB[bid], (a0 - e0) if side == "L" else (s0 - a1))
                break
    if "L" in near and "R" in near:
        (_, l0, l1), _ = near["L"]
        (_, r0, r1), _ = near["R"]
        lo, hi = min(l0, l1, r0, r1), max(l0, l1, r0, r1)
    else:
        (_, g0, g1), d = next(iter(near.values()))
        lo, hi = g0 - 2 * d, g1 + 2 * d
    B = dict(read_fasta(gB))[bseq]
    lo, hi = max(0, lo), min(len(B), hi)
    inside = "-"
    if locB != "none":
        sB, b0, b1 = parse_locus(locB)
        inside = "yes" if (sB == bseq and lo <= b0 and b1 <= hi) else "no"
    best, e, ok, n = search(B, lo, hi, bseq, unitB, hmm, cm, out)
    print(f"{len(on)}\t{bseq}\t{lo}-{hi}\t{(hi - lo) / 1000:.1f}\t{inside}\t{best}\t{e}\t{ok}\t{n}")


if __name__ == "__main__":
    main()
