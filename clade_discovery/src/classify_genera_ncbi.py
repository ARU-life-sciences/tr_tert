#!/usr/bin/env python3
"""Classifies a list of genus names into NCBI's informal higher-clade
label ("blast_name"/"group_name" - e.g. "spiders", "beetles",
"centipedes", "moths & butterflies") via the NCBI Datasets taxonomy API.

This is how clade_discovery picks which DToL genera belong to a target
clade (e.g. Araneae, Myriapoda, Coleoptera) without hand-classifying
hundreds of genus names from memory - NCBI's group_name field is a
reliable, cheap proxy for "which informal higher group is this genus in."

Caches every lookup to a local TSV (default:
clade_discovery/inputs/genus_ncbi_classification.tsv) so re-runs (e.g.
adding a new candidate clade later) don't repeat API calls for genera
already classified.

Usage:
  classify_genera_ncbi.py <genus_list_file> [cache_tsv]

  genus_list_file: one genus name per line (e.g. from `sed -E
  's/_.*//'` on a DToL species-directory listing).
"""
import sys
import time
import urllib.request
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{}/name_report"


def load_cache(path):
    cache = {}
    if path.exists():
        with open(path) as fh:
            next(fh, None)  # header
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 2:
                    cache[parts[0]] = parts[1]
    return cache


def lookup(genus):
    url = API.format(genus)
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.load(resp)
        reports = data.get("reports")
        if not reports:
            return "NOT_FOUND"
        tax = reports[0]["taxonomy"]
        if tax.get("rank") != "GENUS":
            return f"NOT_GENUS_RANK:{tax.get('rank', '?')}"
        return tax.get("group_name", "UNKNOWN")
    except Exception as e:
        return f"ERROR:{e}"


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    genus_list_file = sys.argv[1]
    cache_path = Path(sys.argv[2]) if len(sys.argv) > 2 else (
        PROJECT_ROOT / "inputs" / "genus_ncbi_classification.tsv")

    with open(genus_list_file) as fh:
        genera = [line.strip() for line in fh if line.strip()]
    genera = sorted(set(genera))

    cache = load_cache(cache_path)
    todo = [g for g in genera if g not in cache]
    print(f"[classify_genera_ncbi] {len(genera)} genera, "
          f"{len(cache)} already cached, {len(todo)} to look up",
          file=sys.stderr)

    for i, g in enumerate(todo, 1):
        cache[g] = lookup(g)
        if i % 25 == 0 or i == len(todo):
            print(f"  ({i}/{len(todo)}) {g}: {cache[g]}", file=sys.stderr)
        time.sleep(0.34)  # stay comfortably under NCBI's ~3 req/s guidance

    with open(cache_path, "w") as fh:
        fh.write("genus\tncbi_group_name\n")
        for g in sorted(cache):
            fh.write(f"{g}\t{cache[g]}\n")
    print(f"[classify_genera_ncbi] wrote {cache_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
