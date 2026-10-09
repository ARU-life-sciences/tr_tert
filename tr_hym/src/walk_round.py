#!/usr/bin/env python3
"""One round of model walking for the Hymenoptera TR search.

Seeds: the published TRs (tr_hym/models/*.fa) plus every DToL TR accepted
in an earlier round. A DToL TR is accepted if its best hit has E <= 1e-10
and its template supports the species' own tidk call (template_check.tsv
agrees = yes): the telomere independently validates it, so a wrong hit is
not fed back into the models.

The new round directory gets:
  models/<group>.fa   seeds grouped by family (published families keep
                      their groups; DToL TRs join their family's group, or
                      a new one at >= 3 TRs; the rest are pooled)
  species.tsv         the Hymenoptera with no accepted TR yet
and is then run like tr_hym (build_models.bash, worker_search.bash,
template_check.py, with CLADE_DIR set to it).

Usage: walk_round.py <new_round_dir> <earlier round dir>...
  e.g. walk_round.py ../tr_hym_r2 ../tr_hym
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

MAX_EVALUE = 1e-10
FLANK = 300
COMP = str.maketrans("ACGTN", "TGCAN")


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def core(header, seq):
    """The hit itself from a hits.fa record (seq:a-b|model|E, extracted
    with FLANK either side, in the TR's orientation)."""
    loc = header.split("|")[0]
    a, b = (int(x) for x in loc.rsplit(":", 1)[1].split("-"))
    if a < b:
        off, n = a - max(1, a - FLANK), b - a + 1
    else:
        off, n = FLANK, a - b + 1
    return seq[off:off + n].upper()


def accepted(round_dir):
    """species -> (family, TR sequence) accepted in this round."""
    out = {}
    tc = round_dir / "outputs" / "template_check.tsv"
    if not tc.exists():
        return out
    for r in csv.DictReader(open(tc), delimiter="\t"):
        if r["agrees"] != "yes" or not r["tr_evalue"] or float(r["tr_evalue"]) > MAX_EVALUE:
            continue
        best = None
        for h, s in read_fasta(round_dir / "outputs" / f"{r['species']}.hits.fa"):
            e = float(h.split("|")[2])
            if best is None or e < best[0]:
                best = (e, h, s)
        if best:
            out[r["species"]] = (r["family"], core(best[1], best[2]))
    return out


def main():
    new = Path(sys.argv[1]).resolve()
    earlier = [Path(d).resolve() for d in sys.argv[2:]]
    base = earlier[0]  # tr_hym: published models and the full species list
    groups = defaultdict(list)
    for fa in sorted((base / "models").glob("*.fa")):
        if fa.name.startswith("all_") or fa.name.endswith(".aln.fa"):
            continue
        for h, s in read_fasta(fa):
            groups[fa.stem].append((h, s))
    acc = {}
    for d in earlier:
        acc.update(accepted(d))
    dtol = defaultdict(list)
    for sp, (fam, s) in sorted(acc.items()):
        dtol[fam].append((f"dtol|{sp}|{fam}", s))
    for fam, items in dtol.items():
        if fam in groups or len(items) >= 3:
            groups[fam] += items
        else:
            groups["other_families"] += items
    (new / "models").mkdir(parents=True, exist_ok=True)
    (new / "outputs").mkdir(exist_ok=True)
    (new / "logs").mkdir(exist_ok=True)
    for g, items in groups.items():
        with open(new / "models" / f"{g}.fa", "w") as fh:
            for h, s in items:
                fh.write(f">{h}\n{s}\n")
    species = list(csv.DictReader(open(base / "species.tsv"), delimiter="\t"))
    todo = [r for r in species if r["species"] not in acc]
    with open(new / "species.tsv", "w") as fh:
        print("\t".join(species[0].keys()), file=fh)
        for r in todo:
            print("\t".join(r.values()), file=fh)
    print(f"[walk_round] {len(acc)} accepted DToL TRs as seeds "
          f"({', '.join(f'{f} {len(v)}' for f, v in sorted(dtol.items(), key=lambda x: -len(x[1])))}); "
          f"{len(todo)} species left to search -> {new}")


if __name__ == "__main__":
    main()
