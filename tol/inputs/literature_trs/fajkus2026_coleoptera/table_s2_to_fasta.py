#!/usr/bin/env python3
"""Beetle TRs from Fajkus et al. 2026 GBE 18:evag069 Table S2 (best
Infernal hit per species, column 'TR hit sequence') -> ../coleoptera_tr.fa
(>Genus_species|family|telomere motif), and the paper's per-species TERT
(tblastn) and telomere calls -> table_s2.tsv."""
import openpyxl
from pathlib import Path
HERE = Path(__file__).resolve().parent
rows = list(openpyxl.load_workbook(HERE / "Tables S1-S2" / "Table_S2.xlsx", read_only=True)
            .worksheets[0].iter_rows(values_only=True))
n = 0
with open(HERE.parent / "coleoptera_tr.fa", "w") as fa, open(HERE / "table_s2.tsv", "w") as tsv:
    print("species\tfamily\tGCA\ttelomere_t2t\ttelomere_trfi\ttert\ttr_accession\ttr_from\ttr_to\ttr_strand\ttr_copies", file=tsv)
    for r in rows[6:]:
        if not r[3]:
            continue
        sp = str(r[3]).strip().replace(" ", "_")
        cells = ["" if x is None else str(x).strip() for x in r[:15]]
        print(sp, cells[1], cells[4], cells[6], cells[7], cells[8], cells[9], cells[10], cells[11], cells[12], cells[14], sep="\t", file=tsv)
        seq = "".join(c for c in cells[13].upper() if c in "ACGTU").replace("U", "T")
        if seq:
            fa.write(f">{sp}|{cells[1]}|{cells[7]}\n{seq}\n")
            n += 1
print(f"[table_s2_to_fasta] {n} TR sequences -> coleoptera_tr.fa")
