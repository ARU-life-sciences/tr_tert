#!/usr/bin/env python3
"""One TERT call per species from all the searches, and a table per lineage.

  1. miniprot, three rounds with a growing panel (outputs/tert_round{2,3}_summary.tsv;
     round 1 is superseded by round 2, which searched every species): the best
     status over rounds, present > disrupted > partial > absent
  2. for species still absent: the profile-HMM search of six-frame ORFs
     (outputs/tert_hmm/), where TERT is called if an ORF hits Telomerase_RBD
     at i-Evalue <= HMM_RBD_EVALUE. Calibration: no Diptera (0/6, no TERT
     expected) hit at that level; 6/8 species with miniprot TERTs do.

Calls:
  found        miniprot "present": >= 60% of a reference, intact
  found_fs     miniprot "disrupted": >= 60% of a reference, with 1-3
               frameshifts or stops (mostly alignment or assembly errors:
               35 controls, with working telomeres, are in this class)
  found_hmm    RBD-domain hit in the HMM search only
  partial      miniprot alignment to < 60% of a reference, nothing better
  not_found    none of the above

Outputs: outputs/tert_calls.tsv (per species), outputs/tert_lineages.tsv
"""
import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HMM_RBD_EVALUE = 1e-4
RANK = {"present": 3, "disrupted": 2, "partial": 1, "absent": 0}
CALL = {"present": "found", "disrupted": "found_fs", "partial": "partial", "absent": "not_found"}


def best_rbd(path):
    best = None
    if path.exists():
        for line in open(path):
            if not line.startswith("#"):
                e = float(line.split()[12])
                best = e if best is None else min(best, e)
    return best


def main():
    species = list(csv.DictReader(open(ROOT / "inputs" / "species.tsv"), delimiter="\t"))
    rounds = {}
    for rnd in ("tert_round2", "tert_round3"):
        for r in csv.DictReader(open(ROOT / "outputs" / f"{rnd}_summary.tsv"), delimiter="\t"):
            if r["status"] == "not_run":
                continue
            prev = rounds.get(r["species"])
            if prev is None or RANK[r["status"]] > RANK[prev[0]["status"]]:
                rounds[r["species"]] = (r, rnd)
    cols = ["species", "role", "lineage", "unit", "order", "family", "call", "source",
            "best_query", "best_cover", "frameshifts_stops", "hmm_rbd_evalue"]
    out = []
    for s in species:
        r, rnd = rounds[s["species"]]
        row = {k: s[k] for k in ("species", "role", "lineage", "unit", "order", "family")}
        row.update(call=CALL[r["status"]], source=f"miniprot {rnd[-6:]}",
                   best_query=r["best_query"], best_cover=r["best_cover"],
                   frameshifts_stops=(int(r["frameshifts"] or 0) + int(r["stops"] or 0))
                   if r["best_query"] else "")
        if row["call"] == "not_found":
            e = best_rbd(ROOT / "outputs" / "tert_hmm" / f"{s['species']}.rbd.domtbl")
            row["hmm_rbd_evalue"] = "" if e is None else f"{e:.1e}"
            row["source"] = "miniprot + HMM"
            if e is not None and e <= HMM_RBD_EVALUE:
                row["call"] = "found_hmm"
        out.append(row)
    with open(ROOT / "outputs" / "tert_calls.tsv", "w") as fh:
        print("\t".join(cols), file=fh)
        for r in out:
            print("\t".join(str(r.get(c, "")) for c in cols), file=fh)

    by = defaultdict(lambda: defaultdict(Counter))
    for r in out:
        by[r["lineage"]][r["role"]][r["call"]] += 1
    found = ("found", "found_fs", "found_hmm")
    with open(ROOT / "outputs" / "tert_lineages.tsv", "w") as fh:
        print("lineage\tshift_species\tshift_tert_found\tcontrol_species\tcontrol_tert_found"
              "\ttest_species\ttest_tert_found", file=fh)
        for lin, roles in sorted(by.items(), key=lambda x: -sum(x[1]["shift"].values())):
            cells = [lin]
            for role in ("shift", "control", "test"):
                c = roles[role]
                cells += [sum(c.values()), sum(c[k] for k in found)]
            print("\t".join(map(str, cells)), file=fh)
    print(f"[combine_tert] {len(out)} species: {dict(Counter(r['call'] for r in out))}")


if __name__ == "__main__":
    main()
