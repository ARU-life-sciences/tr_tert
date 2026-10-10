#!/usr/bin/env python3
"""One TERT call per species from the miniprot + Telomerase_RBD search
(src/run_tert.bash).

Alignments of the reference panel are grouped into loci (same sequence and
strand, overlapping). A locus is TERT if the Telomerase_RBD domain (Pfam
PF12009) is found in any of its translations at i-Evalue <= RBD_EVALUE;
other loci are RT-like (e.g. retrotransposon reverse transcriptases, which
share TERT's RT domain but have no RBD) and are counted, not called.

For each TERT locus, the representative alignment is the one with the best
miniprot score (AS), and coverage is that of its reference protein.

Species status:
  present   a TERT locus whose best alignment covers >= MIN_COVER of the
            reference, with no frameshifts or in-frame stops
  disrupted the best TERT locus covers >= MIN_COVER but has frameshifts or
            stops (a pseudogene, or an assembly error)
  partial   TERT loci, but none covering >= MIN_COVER (fragmented in the
            assembly, or a divergent gene the panel aligns to only in part)
  absent    no TERT locus

Outputs (prefix OUT_PREFIX, default outputs/tert_):
  outputs/tert_summary.tsv   one row per species in inputs/species.tsv
  outputs/tert_loci.tsv      every TERT locus
  outputs/tert_proteins.faa  the representative translation of each
                             species' best TERT locus
"""
import csv
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
# a search round's outputs (src/run_tert.bash OUT_DIR), and the prefix for
# this script's: TERT_DIR=outputs/tert_round2 OUT_PREFIX=outputs/tert_round2_
TERT_DIR = ROOT / os.environ.get("TERT_DIR", "outputs/tert")
OUT_PREFIX = os.environ.get("OUT_PREFIX", "outputs/tert_")
RBD_EVALUE = 1e-5
MIN_COVER = 0.6


def read_paf(path):
    """Alignments, each with its translation (the ##STA line after it)."""
    alns, last = [], None
    for line in open(path):
        if line.startswith("##STA"):
            if last is not None:
                last["protein"] = line.rstrip("\n").split("\t")[1]
            last = None
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        tags = dict(t.split(":", 2)[::2] for t in f[12:] if t.count(":") >= 2)
        last = dict(
            query=f[0].split("|")[1] if "|" in f[0] else f[0],
            qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4],
            seq=f[5], ts=int(f[7]), te=int(f[8]),
            score=int(tags.get("AS", 0)), positives=int(tags.get("np", 0)),
            frameshifts=int(tags.get("fs", 0)), stops=int(tags.get("st", 0)),
            protein="",
        )
        last["name"] = f"{f[5]}:{f[7]}-{f[8]}:{f[4]}|{f[0]}"
        alns.append(last)
    return alns


def read_rbd(path):
    """Translation names with a Telomerase_RBD domain at i-Evalue <= RBD_EVALUE."""
    hits = set()
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        if float(f[12]) <= RBD_EVALUE:
            hits.add(f[0])
    return hits


def loci(alns):
    """Group alignments by sequence and strand into overlapping intervals."""
    groups = defaultdict(list)
    for a in alns:
        groups[(a["seq"], a["strand"])].append(a)
    out = []
    for (seq, strand), al in groups.items():
        al.sort(key=lambda a: a["ts"])
        cur, end = [], -1
        for a in al:
            if cur and a["ts"] > end:
                out.append(cur)
                cur, end = [], -1
            cur.append(a)
            end = max(end, a["te"])
        if cur:
            out.append(cur)
    return out


def main():
    species = list(csv.DictReader(open(ROOT / "inputs" / "species.tsv"), delimiter="\t"))
    summary_cols = ["species", "role", "lineage", "unit", "order", "family", "status",
                    "tert_loci", "rt_like_loci", "best_query", "best_cover",
                    "best_positives", "frameshifts", "stops", "locus"]
    locus_cols = ["species", "locus", "best_query", "cover", "positives",
                  "frameshifts", "stops", "protein_length", "queries_aligned"]
    summary, locus_rows, proteins = [], [], []
    for s in species:
        sp = s["species"]
        paf, dom = TERT_DIR / f"{sp}.miniprot.paf", TERT_DIR / f"{sp}.rbd.domtbl"
        row = {k: s[k] for k in ("species", "role", "lineage", "unit", "order", "family")}
        if not dom.exists():
            row["status"] = "not_run"
            summary.append(row)
            continue
        alns, rbd = read_paf(paf), read_rbd(dom)
        tert, rt_like = [], 0
        for group in loci(alns):
            if any(a["name"] in rbd for a in group):
                best = max(group, key=lambda a: a["score"])
                tert.append((best, len({a["query"] for a in group})))
            else:
                rt_like += 1
        for best, nq in tert:
            locus_rows.append(dict(
                species=sp, locus=f"{best['seq']}:{best['ts']}-{best['te']}:{best['strand']}",
                best_query=best["query"], cover=f"{(best['qe'] - best['qs']) / best['qlen']:.2f}",
                positives=best["positives"], frameshifts=best["frameshifts"],
                stops=best["stops"], protein_length=len(best["protein"]), queries_aligned=nq))
        row.update(tert_loci=len(tert), rt_like_loci=rt_like)
        if not tert:
            row["status"] = "absent"
        else:
            best, _ = max(tert, key=lambda t: t[0]["score"])
            cover = (best["qe"] - best["qs"]) / best["qlen"]
            if cover < MIN_COVER:
                row["status"] = "partial"
            elif best["frameshifts"] or best["stops"]:
                row["status"] = "disrupted"
            else:
                row["status"] = "present"
            row.update(best_query=best["query"], best_cover=f"{cover:.2f}",
                       best_positives=best["positives"], frameshifts=best["frameshifts"],
                       stops=best["stops"],
                       locus=f"{best['seq']}:{best['ts']}-{best['te']}:{best['strand']}")
            proteins.append((sp, row["locus"], best["query"], best["protein"]))
        summary.append(row)

    with open(ROOT / f"{OUT_PREFIX}summary.tsv", "w") as fh:
        print("\t".join(summary_cols), file=fh)
        for r in summary:
            print("\t".join(str(r.get(c, "")) for c in summary_cols), file=fh)
    with open(ROOT / f"{OUT_PREFIX}loci.tsv", "w") as fh:
        print("\t".join(locus_cols), file=fh)
        for r in locus_rows:
            print("\t".join(str(r[c]) for c in locus_cols), file=fh)
    with open(ROOT / f"{OUT_PREFIX}proteins.faa", "w") as fh:
        for sp, locus, query, prot in proteins:
            print(f">{sp} {locus} best_query={query}\n{prot}", file=fh)

    counts = defaultdict(lambda: defaultdict(int))
    for r in summary:
        counts[r["role"] if r["role"] != "test" else r["lineage"]][r["status"]] += 1
    for k, v in sorted(counts.items()):
        print(f"[summarise_tert] {k}: {dict(v)}")


if __name__ == "__main__":
    main()
