#!/usr/bin/env python3
"""Species pairs for template co-variation (src/covariation.py): congeners
with different strongly called telomeric units.

Units must look telomeric: <= MAX_UNIT nt and containing CC (the C-rich
strand), which leaves out satellite-like calls (AAAT, ACAGAG,
AAAAGGCAAG, 30-mers), and pairs whose units are the same primitive
unit (a multimer call) are not changes. In each pair, species A (whose template-like
windows are the candidates) has the longer unit, as short units (AACCT:
7 nt windows) give far more chance matches; ties go to the rarer unit in
the genus. Each A is paired with up to MAX_B congeners with another unit.

Inputs: the database calls (telo_calls.tsv) and the ToL assembly list.
Output: pairs.tsv  genus  family  order  A  unit_A  fasta_A  B  unit_B  fasta_B
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path
import sys

sys.path.insert(0, "/lustre/scratch122/tol/teams/blaxter/users/mb39/a-telomeric-repeat-database/src")
from import_tidk_v03 import primitive_unit  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CALLS = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/a-telomeric-repeat-database/data/tidk_v0.3/telo_calls.tsv")
ASSEMBLIES = ROOT.parent.parent / "inputs" / "tol_latest_assemblies.txt"
MAX_UNIT = 16
MAX_B = 3

path = {}
for line in open(ASSEMBLIES):
    sp, p = line.rstrip("\n").split("\t")
    path[p.split("/curated/")[1].split("/")[0]] = (sp, p)

genera = defaultdict(list)
for r in csv.DictReader(open(CALLS), delimiter="\t"):
    u = r["telomere_unit"]
    if r["support"] != "strong" or not r["genus"] or r["assembly"] not in path:
        continue
    if len(u) > MAX_UNIT or "CC" not in u:
        continue
    sp, p = path[r["assembly"]]
    genera[(r["genus"], r["family"], r["order"])].append((sp, u, p))

pairs = []
for (genus, family, order), members in sorted(genera.items()):
    counts = Counter(u for _, u, _ in members)
    if len(counts) < 2:
        continue
    for a, ua, pa in members:
        # a multimer of the same unit (AAACAAACCCAAACCC) is not a change
        others = [m for m in members if primitive_unit(m[1]) != primitive_unit(ua)]
        # A supplies the candidates if its unit is longer, or as long and rarer
        others = [m for m in others
                  if (len(ua), -counts[ua], ua) > (len(m[1]), -counts[m[1]], m[1])]
        for b, ub, pb in sorted(others)[:MAX_B]:
            pairs.append((genus, family, order, a, ua, pa, b, ub, pb))

with open(ROOT / "pairs.tsv", "w") as fh:
    print("genus\tfamily\torder\tA\tunit_A\tfasta_A\tB\tunit_B\tfasta_B", file=fh)
    for p in pairs:
        print("\t".join(p), file=fh)
print(f"[make_pairs] {len(pairs)} pairs in {len({p[0] for p in pairs})} genera")
for p in pairs:
    print(p[0], p[3], p[4], "vs", p[6], p[7])
