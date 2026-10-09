#!/usr/bin/env python3
"""Pairs for the TR-synteny survey: is a TR's genomic position conserved
between two species, by how closely they are related?

TR loci: every TR whose template supports its species' tidk call -
Lepidoptera (../tr_lep, E <= 1e-10) and Hymenoptera (model-walking rounds
../tr_hym..r6, E <= 1e-10, plus the CM/synteny-confirmed sawflies).
Taxonomy: NCBI lineage (taxidlineage.dmp + nodes.dmp ranks), via each
species' taxon_id in the database calls.

A pair's level is the lowest rank the two species share: genus, subfamily,
family, superfamily, order (they differ at every rank below it). Up to
N_PER_LEVEL pairs per clade and level, sampled at random (seed 1), each
species at most twice as A.

Output: survey_pairs.tsv  clade level A A_TR fasta_A B B_TR fasta_B unit_B busco_lineage
"""
import csv
import random
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent      # shifts/denovo
SHIFTS = ROOT.parent
TOL = SHIFTS.parent
TAX = Path("/data/tol/resources/taxonomy/2026-10-03/new_taxdump")
CALLS = Path("/lustre/scratch122/tol/teams/blaxter/users/mb39/a-telomeric-repeat-database/data/tidk_v0.3/telo_calls.tsv")
LEVELS = ["genus", "subfamily", "family", "superfamily", "order"]
N_PER_LEVEL = 25
MAX_EVALUE = 1e-10
EXTRA_HYM = {  # confirmed by CM + synteny + template (../tr_cm, README)
    "Allantus_rufocinctus": "SUPER_4:15750949-15750823",
    "Abia_aenea": "SUPER_7:6670205-6670081",
    "Athalia_cordata": "SUPER_1:40834657-40834534",
    "Athalia_circularis": "SUPER_4:5800764-5800695",
}


def loci():
    """species -> (clade, TR locus, unit, fasta)."""
    out = {}
    fasta = {}
    for d in ("tr_hym", "tr_lep"):
        for r in csv.DictReader(open(SHIFTS / d / "species.tsv"), delimiter="\t"):
            fasta[r["species"]] = r["fasta"]
    for clade, dirs in (("Lepidoptera", ["tr_lep"]),
                        ("Hymenoptera", ["tr_hym"] + [f"tr_hym_r{i}" for i in range(2, 7)])):
        for d in dirs:
            for r in csv.DictReader(open(SHIFTS / d / "outputs" / "template_check.tsv"), delimiter="\t"):
                if r["agrees"] == "yes" and r["tr_evalue"] and float(r["tr_evalue"]) <= MAX_EVALUE:
                    out.setdefault(r["species"], (clade, r["tr_locus"], r["tidk_unit"], fasta[r["species"]]))
    for sp, loc in EXTRA_HYM.items():
        out.setdefault(sp, ("Hymenoptera", loc, "AACCT", fasta[sp]))
    return out


def lineages(taxids):
    rank, name = {}, {}
    with open(TAX / "nodes.dmp") as fh:
        for line in fh:
            f = line.split("\t|\t")
            rank[f[0]] = f[2]
    lin = {}
    with open(TAX / "taxidlineage.dmp") as fh:
        for line in fh:
            t, rest = line.split("\t|\t", 1)
            if t in taxids:
                lin[t] = rest.strip("\t|\n").split() + [t]
    wanted = {x for l in lin.values() for x in l if rank.get(x) in LEVELS}
    with open(TAX / "names.dmp") as fh:
        for line in fh:
            f = line.split("\t|\t")
            if f[0] in wanted and f[3].startswith("scientific name"):
                name[f[0]] = f[1]
    return {t: {rank[x]: name.get(x, x) for x in l if rank.get(x) in LEVELS} for t, l in lin.items()}


def main():
    L = loci()
    tolid = {fa.split("/curated/")[1].split("/")[0]: sp for sp, (_, _, _, fa) in L.items()}
    taxid = {}
    for r in csv.DictReader(open(CALLS), delimiter="\t"):
        if r["assembly"] in tolid and r["taxon_id"]:
            taxid[tolid[r["assembly"]]] = r["taxon_id"]
    lin = lineages(set(taxid.values()))
    tax = {sp: lin.get(t, {}) for sp, t in taxid.items()}
    random.seed(1)
    rows = []
    for clade in ("Lepidoptera", "Hymenoptera"):
        sps = sorted(s for s in L if L[s][0] == clade and tax.get(s))
        by_level = defaultdict(list)
        for i, a in enumerate(sps):
            for b in sps[i + 1:]:
                ta, tb = tax[a], tax[b]
                for lev in LEVELS:
                    if ta.get(lev) and ta.get(lev) == tb.get(lev):
                        lower = LEVELS[:LEVELS.index(lev)]
                        if all(ta.get(x) != tb.get(x) or not ta.get(x) for x in lower):
                            by_level[lev].append((a, b))
                        break
        lineage = "lepidoptera_odb12" if clade == "Lepidoptera" else "hymenoptera_odb12"
        for lev in LEVELS:
            pairs = by_level[lev][:]
            random.shuffle(pairs)
            used = defaultdict(int)
            n = 0
            for a, b in pairs:
                if n == N_PER_LEVEL:
                    break
                if used[a] >= 2:
                    continue
                used[a] += 1
                n += 1
                rows.append([clade, lev, a, L[a][1], L[a][3], b, L[b][1], L[b][3], L[b][2], lineage])
            print(f"[survey_pairs] {clade} {lev}: {len(by_level[lev])} possible, {n} sampled")
    with open(ROOT / "survey_pairs.tsv", "w") as fh:
        print("clade\tlevel\tA\tA_TR\tfasta_A\tB\tB_TR\tfasta_B\tunit_B\tbusco_lineage", file=fh)
        for r in rows:
            print("\t".join(r), file=fh)


if __name__ == "__main__":
    main()
