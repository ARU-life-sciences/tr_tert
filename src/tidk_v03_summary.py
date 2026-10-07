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
  tidk_v03_summary.py [outfile]    # default <EXPLORE_DIR>/summary.tsv

Paths default to the plant run; override with environment variables for
other runs (e.g. clade_discovery):
  ASSEMBLY_LIST  species<TAB>fasta list(s), colon-separated
                 (default inputs/dtol_plant_paths.txt)
  OLD_TIDK_DIR   tidk 0.2.7 results (default outputs/tidk)
  EXPLORE_DIR    tidk 0.3.0 explore results (default outputs/tidk_v0.3)
  ENDS_DIR       tidk ends results (default outputs/tidk_ends)

Columns:
  species  v027_top  v027_top_unit  v03_top  v03_copies  v03_minor_strand
  relation  telomere_unit  telomere_rank  telomere_minor_strand  sequences
  t2t  one_end  no_ends  wrong_strand  not_terminal

v03_minor_strand is for the top explore unit, telomere_minor_strand for
the telomere call (the same unit unless telomere_rank > 1).

relation compares the top repeats of the two versions:
  same       the v0.2.7 top reduces to exactly the 0.3.0 unit (it only lost
             multimer padding)
  compound   the v0.2.7 top contains a full copy of the 0.3.0 unit, but is a
             mixed/compound string (e.g. AAAACCCTAAACCCT for AAACCCT)
  different  neither: the call changed
  none       one of the versions found nothing
"""
import json
import os
import sys
from pathlib import Path

import tidk_io

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "outputs"


def env_path(name, default):
    p = Path(os.environ.get(name, default))
    return p if p.is_absolute() else PROJECT_ROOT / p


ASSEMBLY_LISTS = [Path(l) if Path(l).is_absolute() else PROJECT_ROOT / l
                  for l in os.environ.get("ASSEMBLY_LIST", "inputs/dtol_plant_paths.txt").split(":")]
OLD_TIDK_DIR = env_path("OLD_TIDK_DIR", "outputs/tidk")
EXPLORE_DIR = env_path("EXPLORE_DIR", "outputs/tidk_v0.3")
ENDS_DIR = env_path("ENDS_DIR", "outputs/tidk_ends")

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


def load_explore_rows(species):
    """{unit: (copies, minor strand proportion)} and the units in rank
    order, from the 0.3.0 explore output; None if there isn't one."""
    path = EXPLORE_DIR / f"{species}.explore.tsv"
    if not path.exists():
        return None
    rows, order = {}, []
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            unit, copies, as_unit, as_revcomp = parts[:4]
            minor = None
            if as_unit != "NA" and int(as_unit) + int(as_revcomp) > 0:
                u, r = int(as_unit), int(as_revcomp)
                minor = min(u, r) / (u + r)
            rows[unit] = (int(copies), minor)
            order.append(unit)
    return rows, order


def fmt_minor(minor):
    return "NA" if minor is None else f"{minor:.3f}"


def present_ends(ends):
    return 2 * ends["t2t"] + ends["one_end"]


def telomere_call(species):
    """(explore rank, ends summary) of the candidate unit at the most
    chromosome ends; ties go to the higher explore rank. Falls back to the
    main run's ends for the top explore unit."""
    cands = []
    for path in (ENDS_DIR / "candidates").glob(f"{species}.*.ends.json"):
        rank = path.name[len(species) + 1:].split(".")[0]
        if rank.isdigit():
            cands.append((int(rank), json.load(open(path))))
    if cands:
        return max(cands, key=lambda c: (present_ends(c[1]), -c[0]))
    path = ENDS_DIR / f"{species}.ends.json"
    if path.exists():
        return (1, json.load(open(path)))
    return (None, {})


def main():
    outfile = sys.argv[1] if len(sys.argv) > 1 else EXPLORE_DIR / "summary.tsv"
    cols = ["species", "v027_top", "v027_top_unit", "v03_top", "v03_copies",
            "v03_minor_strand", "relation", "telomere_unit", "telomere_rank",
            "telomere_minor_strand", "sequences", "t2t",
            "one_end", "no_ends", "wrong_strand", "not_terminal"]
    rows = []
    species_list = [line.split("\t")[0].strip() for l in ASSEMBLY_LISTS for line in open(l)]
    for species in species_list:
        if not species:
            continue
        explore = load_explore_rows(species)
        if explore is None:
            continue
        rows_by_unit, order = explore
        top = order[0] if order else ""
        v03 = (top, *rows_by_unit[top]) if top else ("", 0, None)
        _, old = tidk_io.load_explore(OLD_TIDK_DIR / f"{species}.tidk.tsv")
        old_top = old[0][0] if old else ""
        old_unit = primitive_unit(old_top) if old_top else ""
        rank, ends = telomere_call(species)
        minor = v03[2]
        rows.append(dict(
            species=species, v027_top=old_top, v027_top_unit=old_unit,
            v03_top=v03[0], v03_copies=v03[1],
            v03_minor_strand=fmt_minor(minor),
            relation=relation(old_top, v03[0]),
            telomere_unit=ends.get("unit", ""), telomere_rank="" if rank is None else rank,
            telomere_minor_strand=fmt_minor(rows_by_unit.get(ends.get("unit", ""), (0, None))[1]),
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
