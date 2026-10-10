#!/usr/bin/env python3
"""Pick the species for the TERT/TR search: the animal telomeric-repeat
shift lineages, matched controls, and groups that may lack telomerase.

A shift lineage is a non-ancestral telomere call (strong support: >= 25% of
chromosome ends) shared by >= 2 species of one order, outside Diptera
(which lost telomerase; their calls are satellites at the ends). Ancestral units:
AACCT for arthropods, AACCCT for other animals. Plants are left out: their
TR/TERT are in ../../tr_tert.

For each lineage: every shift species, plus up to N_CONTROLS controls with
the ancestral unit (strong), from the shift species' families first, then
the rest of the order, in species order (so the choice is reproducible).

Test groups (TERT may be absent): Diptera (telomerase lost in Drosophila
and relatives), Araneae (no telomere at any chromosome end in 18/18), and
Ephemeroptera (no call in 8/8), up to N_TESTS species each, preferring the
most chromosome ends checked.

Input:  ../../a-telomeric-repeat-database/data/tidk_v0.3/telo_calls.tsv
        ../inputs/tol_latest_assemblies.txt
Output: inputs/species.tsv
        species  fasta  role (shift/control/test)  lineage  unit  support
        order  family
"""
import csv
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
CALLS = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/"
             "a-telomeric-repeat-database/data/tidk_v0.3/telo_calls.tsv")
ASSEMBLIES = HERE.parent / "inputs" / "tol_latest_assemblies.txt"
N_CONTROLS = 5
N_TESTS = 6
TESTS = {"Diptera": "order", "Araneae": "order", "Ephemeroptera": "order"}


def ancestral(r):
    if r["kingdom"] != "Metazoa":
        return None
    return "AACCT" if r["phylum"] == "Arthropoda" else "AACCCT"


def main():
    fasta = {}
    for line in open(ASSEMBLIES):
        sp, path = line.rstrip("\n").split("\t")
        fasta[sp.replace("_", " ")] = (sp, path)
    rows = [r for r in csv.DictReader(open(CALLS), delimiter="\t")]
    # telo_calls uses NCBI names; match to assembly names on the ToLID
    by_assembly = {}
    for sp, path in fasta.values():
        tolid_dir = path.split("/curated/")[1].split("/")[0]
        by_assembly[tolid_dir] = (sp, path)
    for r in rows:
        r["_sp"], r["_fasta"] = by_assembly.get(r["assembly"], ("", ""))
    rows = [r for r in rows if r["_sp"]]

    strong = [r for r in rows if r["support"] == "strong" and ancestral(r)]
    lineages = defaultdict(list)
    for r in strong:
        # Diptera lost telomerase: their "shifts" are satellites at ends
        if r["telomere_unit"] != ancestral(r) and r["order"] and r["order"] != "Diptera":
            lineages[(r["order"], r["telomere_unit"])].append(r)
    lineages = {k: v for k, v in lineages.items() if len(v) >= 2}

    out, seen = [], set()

    def add(r, role, lineage):
        if r["_sp"] in seen:
            return
        seen.add(r["_sp"])
        out.append([r["_sp"], r["_fasta"], role, lineage, r["telomere_unit"],
                    r["support"], r["order"], r["family"]])

    for (order, unit), shift in sorted(lineages.items(), key=lambda x: (-len(x[1]), x[0])):
        lineage = f"{order}:{unit}"
        for r in sorted(shift, key=lambda r: r["_sp"]):
            add(r, "shift", lineage)
        fams = {r["family"] for r in shift}
        pool = [r for r in strong if r["order"] == order and r["telomere_unit"] == ancestral(r)]
        pool.sort(key=lambda r: (r["family"] not in fams, r["_sp"]))
        n = 0
        for r in pool:
            if n == N_CONTROLS:
                break
            if r["_sp"] not in seen:
                add(r, "control", lineage)
                n += 1
    for name, rank in TESTS.items():
        pool = [r for r in rows if r[rank] == name]
        pool.sort(key=lambda r: (-int(r["ends_checked"] or 0), r["_sp"]))
        for r in pool[:N_TESTS]:
            add(r, "test", f"{name}:test")

    with open(HERE / "inputs" / "species.tsv", "w") as fh:
        print("species\tfasta\trole\tlineage\tunit\tsupport\torder\tfamily", file=fh)
        for row in out:
            print("\t".join(row), file=fh)
    roles = defaultdict(int)
    for row in out:
        roles[row[2]] += 1
    print(f"[select_species] {len(out)} species, {len(lineages)} shift lineages: {dict(roles)}")


if __name__ == "__main__":
    main()
