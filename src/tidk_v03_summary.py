#!/usr/bin/env python3
"""Per-species summary of the tidk 0.3.0 rerun (src/run_tidk_v03.bash),
compared with the tidk 0.2.7 results the existing notes are built on.

For each species in the assembly list:
  - the top v0.2.7 repeat, and its primitive canonical unit (rotations,
    reverse complement and exact multimers collapsed, as 0.3.0 does), so
    "the same telomere, minus multimer padding" can be told apart from a
    genuinely different call
  - the top 0.3.0 unit, its copies and the proportion on the minor strand
    (a telomeric repeat reads both ways round, at the two chromosome ends)
  - the telomere call: of the top explore units checked with tidk ends
    (src/run_tidk_ends_candidates.bash), the one at the most chromosome
    ends. The top explore unit isn't always the telomere - a long
    subterminal satellite can have more copies near the ends - so this,
    not v03_top, is the species' telomeric repeat. Falls back to the main
    run's ends (top explore unit) if candidates haven't been run.
  - tidk ends for the telomere call: sequences checked, T2T / one end /
    no ends, and flagged ends

Usage:
  tidk_v03_summary.py [outfile]    # default outputs/tidk_v0.3/summary.tsv

Columns:
  species  v027_top  v027_top_unit  v03_top  v03_copies  v03_minor_strand
  relation  telomere_unit  telomere_rank  sequences  t2t  one_end  no_ends
  wrong_strand  not_terminal

relation compares the top repeats of the two versions:
  same       the v0.2.7 top reduces to exactly the 0.3.0 unit (it only lost
             multimer padding)
  compound   the v0.2.7 top contains a full copy of the 0.3.0 unit, but is a
             mixed/compound string (e.g. AAAACCCTAAACCCT for AAACCCT)
  different  neither: the call changed
  none       one of the versions found nothing
"""
import json
import sys
from pathlib import Path

import tidk_io

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def lex_min(s):
    """Lexicographically smallest rotation of s or its reverse complement."""
    rc = s.translate(COMPLEMENT)[::-1]
    return min(min(x[i:] + x[:i] for i in range(len(x))) for x in (s, rc))


def primitive_unit(s):
    """As tidk 0.3.0's utils::primitive_telomere_unit: canonicalise, reduce
    to the smallest exact period, canonicalise again."""
    canon = lex_min(s)
    n = len(canon)
    for p in range(1, n + 1):
        if n % p == 0 and canon == canon[:p] * (n // p):
            return lex_min(canon[:p])
    return canon


def relation(old_top, new_top):
    if not old_top or not new_top:
        return "none"
    if primitive_unit(old_top) == new_top:
        return "same"
    rc = old_top.translate(COMPLEMENT)[::-1]
    copies = {new_top[i:] + new_top[:i] for i in range(len(new_top))}
    if any(c in s for c in copies for s in (old_top, rc)):
        return "compound"
    return "different"


def load_v03(species):
    """(unit, copies, minor strand proportion) of the top 0.3.0 unit."""
    path = OUTPUT_DIR / "tidk_v0.3" / f"{species}.explore.tsv"
    if not path.exists():
        return None
    with open(path) as fh:
        next(fh, None)
        line = fh.readline().rstrip("\n")
    if not line:
        return ("", 0, None)
    unit, copies, as_unit, as_revcomp = line.split("\t")[:4]
    if as_unit == "NA":
        return (unit, int(copies), None)
    u, r = int(as_unit), int(as_revcomp)
    return (unit, int(copies), min(u, r) / (u + r) if u + r else None)


def present_ends(ends):
    return 2 * ends["t2t"] + ends["one_end"]


def telomere_call(species):
    """(explore rank, ends summary) of the candidate unit at the most
    chromosome ends; ties go to the higher explore rank. Falls back to the
    main run's ends for the top explore unit."""
    cands = []
    for path in (OUTPUT_DIR / "tidk_ends" / "candidates").glob(f"{species}.*.ends.json"):
        rank = path.name[len(species) + 1:].split(".")[0]
        if rank.isdigit():
            cands.append((int(rank), json.load(open(path))))
    if cands:
        return max(cands, key=lambda c: (present_ends(c[1]), -c[0]))
    path = OUTPUT_DIR / "tidk_ends" / f"{species}.ends.json"
    if path.exists():
        return (1, json.load(open(path)))
    return (None, {})


def main():
    outfile = sys.argv[1] if len(sys.argv) > 1 else OUTPUT_DIR / "tidk_v0.3" / "summary.tsv"
    cols = ["species", "v027_top", "v027_top_unit", "v03_top", "v03_copies",
            "v03_minor_strand", "relation", "telomere_unit", "telomere_rank", "sequences", "t2t",
            "one_end", "no_ends", "wrong_strand", "not_terminal"]
    rows = []
    for line in open(PROJECT_ROOT / "inputs" / "dtol_plant_paths.txt"):
        species = line.split("\t")[0].strip()
        if not species:
            continue
        v03 = load_v03(species)
        if v03 is None:
            continue
        _, old = tidk_io.load_explore(OUTPUT_DIR / "tidk" / f"{species}.tidk.tsv")
        old_top = old[0][0] if old else ""
        old_unit = primitive_unit(old_top) if old_top else ""
        rank, ends = telomere_call(species)
        minor = v03[2]
        rows.append(dict(
            species=species, v027_top=old_top, v027_top_unit=old_unit,
            v03_top=v03[0], v03_copies=v03[1],
            v03_minor_strand="NA" if minor is None else f"{minor:.3f}",
            relation=relation(old_top, v03[0]),
            telomere_unit=ends.get("unit", ""), telomere_rank="" if rank is None else rank,
            sequences=ends.get("sequences", ""),
            t2t=ends.get("t2t", ""), one_end=ends.get("one_end", ""),
            no_ends=ends.get("no_ends", ""),
            wrong_strand=len(ends.get("wrong_strand", [])) if ends else "",
            not_terminal=len(ends.get("not_terminal", [])) if ends else "",
        ))
    with open(outfile, "w") as out:
        print("\t".join(cols), file=out)
        for r in rows:
            print("\t".join(str(r[c]) for c in cols), file=out)
    print(f"[tidk_v03_summary] {len(rows)} species -> {outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
