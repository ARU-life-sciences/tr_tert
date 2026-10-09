#!/usr/bin/env python3
"""Summary of the TR-synteny survey (src/worker_survey.bash outputs).

Per pair: conserved if B's TR is inside the interval projected from A's TR
by nucleotide or BUSCO anchors; moved if it is in neither but the anchors
are good (>= MIN_NUC nucleotide or >= MIN_BUSCO gene anchors); otherwise
undetermined. Output: outputs/survey_pairs_result.tsv and
outputs/survey_summary.tsv (clade x level).
"""
import csv
import glob
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MIN_NUC, MIN_BUSCO = 10, 3
LEVELS = ["genus", "subfamily", "family", "superfamily", "order"]


def num(x):
    try:
        return int(x)
    except ValueError:
        return 0


rows = []
for f in glob.glob(str(ROOT / "outputs" / "survey" / "*.tsv")):
    for line in open(f):
        clade, level, a, b, na, ni, ba, bi = (line.rstrip("\n").split("\t") + [""] * 8)[:8]
        if ni == "yes" or bi == "yes":
            call = "conserved"
        elif num(na) >= MIN_NUC or num(ba) >= MIN_BUSCO:
            call = "moved"
        else:
            call = "undetermined"
        rows.append(dict(clade=clade, level=level, A=a, B=b, nuc_anchors=na, nuc_inside=ni,
                         busco_anchors=ba, busco_inside=bi, call=call))
cols = list(rows[0].keys())
with open(ROOT / "outputs" / "survey_pairs_result.tsv", "w") as fh:
    print("\t".join(cols), file=fh)
    for r in sorted(rows, key=lambda r: (r["clade"], LEVELS.index(r["level"]), r["A"])):
        print("\t".join(r[c] for c in cols), file=fh)
summ = defaultdict(Counter)
for r in rows:
    summ[(r["clade"], r["level"])][r["call"]] += 1
with open(ROOT / "outputs" / "survey_summary.tsv", "w") as fh:
    print("clade\tlevel\tpairs\tconserved\tmoved\tundetermined\tconserved_of_determined", file=fh)
    for clade in ("Lepidoptera", "Hymenoptera"):
        for lev in LEVELS:
            c = summ[(clade, lev)]
            det = c["conserved"] + c["moved"]
            frac = f"{c['conserved'] / det:.2f}" if det else "-"
            print(clade, lev, sum(c.values()), c["conserved"], c["moved"], c["undetermined"], frac, sep="\t", file=fh)
print(open(ROOT / "outputs" / "survey_summary.tsv").read())
