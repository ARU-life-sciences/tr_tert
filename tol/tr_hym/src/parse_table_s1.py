#!/usr/bin/env python3
"""Fajkus et al. 2023 Table S1 -> table_s1.tsv and per-family TR sets.

The table lists ranks in merged-looking cells: a rank is filled only where
it changes. The first parse carried each rank forward independently, so
rows of a new superfamily with a blank family (Orussoidea, Cephoidea,
Tenthredinoidea, Formicoidea, ...) inherited the previous family - the
sawfly TRs (Athalia, Cephus, Diprion, Tenthredo, ...) were filed under
Ichneumonidae. Here a new superfamily resets the lower ranks, and a blank
family falls back to the superfamily.

Usage: parse_table_s1.py <Table S1.xlsx> <hymenoptera_tr.fa> <out table_s1.tsv> <out models dir>
Writes one <group>.fa per family, however few TRs it has (headers:
accession|species|motif): nhmmer works from a single sequence, and pooling
singletons buried exactly the lineages that need them (Cephus, Orussus).
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import openpyxl


def clean(x):
    return re.sub(r"[\s\xa0]+", " ", str(x or "")).strip()


def main():
    xlsx, fa, out_tsv, out_dir = sys.argv[1:5]
    ws = openpyxl.load_workbook(xlsx, read_only=True).worksheets[0]
    rows = list(ws.iter_rows(values_only=True))[4:]
    ctx = ["", "", "", ""]
    out = []
    for r in rows:
        r = list(r) + [""] * (21 - len(r))
        if clean(r[0]):
            ctx = [clean(r[0]), "", "", ""]
        for i in (1, 2, 3):
            if clean(r[i]):
                ctx[i] = clean(r[i])
                ctx[i + 1:] = [""] * (3 - i)
        if not clean(r[4]):
            continue
        family = ctx[1] or ctx[0]
        out.append([ctx[0], family, ctx[2], ctx[3], clean(r[4]), clean(r[5]), clean(r[8]),
                    re.sub(r"\s", "", clean(r[12])), clean(r[9]), clean(r[10]), clean(r[13])])
    with open(out_tsv, "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["superfamily", "family", "subfamily", "tribe", "species", "assembly", "motif",
                    "tr_accession", "trfi", "tr_template", "type3_promoter"])
        w.writerows(out)
    meta = {r[7]: r for r in out if r[7].startswith("BK")}
    seqs, name = {}, None
    for line in open(fa):
        if line.startswith(">"):
            name = line[1:].split()[0].split(".")[0]
            seqs[name] = ""
        else:
            seqs[name] += line.strip()
    groups = defaultdict(list)
    for acc, s in seqs.items():
        r = meta[acc]
        groups[r[1]].append((f"{acc}|{r[4].replace(' ', '_')}|{r[6]}", s))
    d = Path(out_dir)
    d.mkdir(parents=True, exist_ok=True)
    for g, items in sorted(groups.items()):
        with open(d / f"{g}.fa", "w") as fh:
            for h, s in items:
                fh.write(f">{h}\n{s}\n")
    print("[parse_table_s1] " + ", ".join(f"{g} {len(v)}" for g, v in
                                            sorted(groups.items(), key=lambda x: -len(x[1]))))


if __name__ == "__main__":
    main()
