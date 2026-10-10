#!/usr/bin/env python3
"""Iterative search: add the intact TERTs found in DToL species (status
"present" in outputs/tert_summary.tsv: >= 60% of a reference, no
frameshifts or stops) to the reference panel, as queries for the next
round. Divergent clades with no panel species (e.g. Carabidae, Adephaga)
are then searched with a relative from their own family.

Usage: expand_panel.py <panel.faa> <tert_summary.tsv> <tert_proteins.faa> <out.faa>
"""
import csv
import sys

panel, summary, proteins, out = sys.argv[1:]
present = {r["species"] for r in csv.DictReader(open(summary), delimiter="\t")
           if r["status"] == "present"}
seqs, name = {}, None
for line in open(proteins):
    if line.startswith(">"):
        name = line[1:].split()[0]
    elif name in present:
        seqs[name] = line.strip()
seen = set()
with open(out, "w") as fh:
    fh.write(open(panel).read())
    for sp, prot in sorted(seqs.items()):
        if prot in seen:
            continue
        seen.add(prot)
        fh.write(f">dtol|{sp}|TERT\n{prot}\n")
print(f"[expand_panel] {len(seen)} DToL TERTs added to {panel} -> {out}")
