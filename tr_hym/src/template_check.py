#!/usr/bin/env python3
"""Does each Hymenoptera species' TR template encode the telomeric repeat
tidk found at its chromosome ends?

For each species (species.tsv), the TR locus is the best nhmmer hit
(lowest E-value; hits from different family models at one locus are the
same locus). The template is read from the hit with its 300 bp flanks, on
the TR's sense strand, where it is ~1.2-2 copies of the telomere's C-rich
strand: e.g. Tetragonula carbonaria (TTTAGGTTAGG) has AAACCTAACCTAAACCT.

For a unit u, the template match is the longest stretch of the locus that
exactly matches the tandem of u in some phase. A unit is supported if
that stretch is >= len(u) + MIN_EXTRA nt: one copy and a little more, as
templates are often ~1.2 units (Bombus: 13 nt for AACCCCAACCT).
agrees = yes if the species' tidk call is supported. If not, the
template's unit is the supported known telomeric unit
(../../inputs/known_telomere_units.tsv) with the longest match, ties to
the shorter unit, and agrees = no; no_template if none is supported.

Output: outputs/template_check.tsv
  species family tidk_unit tidk_support tr_locus tr_model tr_evalue
  tidk_match template_unit template_match agrees
"""
import csv
import os
import re
from pathlib import Path

ROOT = Path(os.environ.get("CLADE_DIR", Path(__file__).resolve().parent.parent))
KNOWN = ROOT.parent.parent / "inputs" / "known_telomere_units.tsv"
MIN_EXTRA = 2
MAX_EVALUE = 1e-10
# PICK_TEMPLATE=1: the TR is the best hit with a supporting template, not
# simply the best hit (see main)
PICK_TEMPLATE = os.environ.get("PICK_TEMPLATE", "0") == "1"


def longest_tandem_match(seq, unit):
    """Longest substring of seq equal to a substring of unit repeated."""
    n = len(unit)
    best = 0
    # from each start, extend while seq follows the tandem from each phase
    for start in range(len(seq)):
        for phase in range(n):
            if seq[start] != unit[phase]:
                continue
            k = 0
            while start + k < len(seq) and seq[start + k] == unit[(phase + k) % n]:
                k += 1
            best = max(best, k)
    return best


def read_hits(path):
    hits = []
    name = None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].strip()
            hits.append([name, ""])
        elif name:
            hits[-1][1] += line.strip().upper()
    out = []
    for name, seq in hits:
        locus, model, evalue = name.split("|")
        out.append((float(evalue), locus, model, seq))
    return sorted(out)


def main():
    known = [r["unit"] for r in csv.DictReader(open(KNOWN), delimiter="\t")]
    rows = []
    for s in csv.DictReader(open(ROOT / "species.tsv"), delimiter="\t"):
        sp = s["species"]
        row = dict(species=sp, family=s["family"], tidk_unit=s["unit"], tidk_support=s["support"])
        path = ROOT / "outputs" / f"{sp}.hits.fa"
        if not path.exists():
            row["tr_locus"] = "not_run"
            rows.append(row)
            continue
        hits = read_hits(path)
        if PICK_TEMPLATE and s["unit"]:
            # genomes can hold several TR-like copies, most with a decayed
            # template (pimplines: 3-7 copies, one with a template): take
            # the best hit whose template supports the tidk call, whatever
            # its E (still reported; validated only at E <= MAX_EVALUE)
            with_t = [h for h in hits if longest_tandem_match(h[3], s["unit"]) >= len(s["unit"]) + MIN_EXTRA]
            hits = with_t or [h for h in hits if h[0] <= MAX_EVALUE]
        else:
            hits = [h for h in hits if h[0] <= MAX_EVALUE]
        if not hits:
            row["tr_locus"] = "none"
            rows.append(row)
            continue
        evalue, locus, model, seq = hits[0]
        row["tr_copies"] = len({h[1].rsplit(":", 1)[0] + ":" + str(int(h[1].rsplit(":", 1)[1].split("-")[0]) // 2000)
                                for h in read_hits(path)})
        row.update(tr_locus=locus, tr_model=model, tr_evalue=f"{evalue:.1e}")
        units = list(dict.fromkeys(([s["unit"]] if s["unit"] else []) + known))
        matches = {u: longest_tandem_match(seq, u) for u in units}
        if s["unit"]:
            row["tidk_match"] = matches[s["unit"]]
        if s["unit"] and matches[s["unit"]] >= len(s["unit"]) + MIN_EXTRA:
            row.update(template_unit=s["unit"], template_match=matches[s["unit"]], agrees="yes")
        else:
            supported = [(m, -len(u), u) for u, m in matches.items() if m >= len(u) + MIN_EXTRA]
            if supported:
                _, _, u = max(supported)
                row.update(template_unit=u, template_match=matches[u], agrees="no")
            else:
                row["agrees"] = "no_template"
        rows.append(row)
    cols = ["species", "family", "tidk_unit", "tidk_support", "tr_locus", "tr_model",
            "tr_evalue", "tidk_match", "template_unit", "template_match", "agrees", "tr_copies"]
    with open(ROOT / "outputs" / "template_check.tsv", "w") as fh:
        print("\t".join(cols), file=fh)
        for r in rows:
            print("\t".join(str(r.get(c, "")) for c in cols), file=fh)
    print(f"[template_check] {len(rows)} species -> outputs/template_check.tsv")


if __name__ == "__main__":
    main()
