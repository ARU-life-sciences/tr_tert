#!/usr/bin/env python3
"""Loci from the sensitive TERT search (src/run_tert_hmm.bash).

Each ORF is one exon or less, so the TERT profile's hits are grouped into
loci: ORFs on the same sequence and strand within MAX_GAP of each other.
A locus is scored by the sum of its ORFs' best domain bit scores, and by
how much of the TERT profile its hits cover (union of HMM coordinates). It
has an RBD if an ORF in it hits Telomerase_RBD at i-Evalue <= RBD_EVALUE.

Retrotransposon reverse transcriptases also hit the profile (TERT's RT
domain is homologous), but only over the RT region and without an RBD:
Diptera, which have no TERT, give the background.

Output: outputs/tert_hmm_summary.tsv, one row per species: the best locus
(by score) and the best locus with an RBD.
"""
import csv
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
HMM_DIR = ROOT / "outputs" / "tert_hmm"
MAX_GAP = 50_000
RBD_EVALUE = 1e-5
COORDS = re.compile(r"source=(\S+) coords=(\d+)\.\.(\d+)")


def read_domtbl(path, max_evalue=None):
    """ORF name -> (sequence, start, end, strand, best bit score,
    [(hmm_from, hmm_to)])."""
    orfs = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split(None, 22)
        ievalue, bits = float(f[12]), float(f[13])
        if max_evalue is not None and ievalue > max_evalue:
            continue
        m = COORDS.search(f[22])
        if not m:
            continue
        seq, a, b = m.group(1), int(m.group(2)), int(m.group(3))
        strand = "+" if a <= b else "-"
        o = orfs.setdefault(f[0], [seq, min(a, b), max(a, b), strand, 0.0, []])
        o[4] = max(o[4], bits)
        o[5].append((int(f[15]), int(f[16])))
    return orfs


def covered(intervals):
    total, end = 0, 0
    for a, b in sorted(intervals):
        if b > end:
            total += b - max(a, end + 1) + 1
            end = b
    return total


def loci(orfs):
    groups = defaultdict(list)
    for name, o in orfs.items():
        groups[(o[0], o[3])].append((name, o))
    out = []
    for key, items in groups.items():
        items.sort(key=lambda x: x[1][1])
        cur = []
        for item in items:
            if cur and item[1][1] - cur[-1][1][2] > MAX_GAP:
                out.append(cur)
                cur = []
            cur.append(item)
        out.append(cur)
    return out


def main():
    hmm_len = int(next(l.split()[1] for l in open(ROOT / "inputs" / "tert_round3.hmm")
                       if l.startswith("LENG")))
    species = list(csv.DictReader(open(ROOT / "inputs" / "species_hmm.tsv"), delimiter="\t"))
    cols = ["species", "role", "lineage", "order", "family", "best_score", "best_cover",
            "best_orfs", "best_rbd", "best_locus", "rbd_score", "rbd_cover", "rbd_locus"]
    rows = []
    for s in species:
        sp = s["species"]
        row = {k: s[k] for k in ("species", "role", "lineage", "order", "family")}
        tert_path = HMM_DIR / f"{sp}.tert.domtbl"
        if not tert_path.exists():
            rows.append(row)
            continue
        tert = read_domtbl(tert_path)
        rbd = set(read_domtbl(HMM_DIR / f"{sp}.rbd.domtbl", RBD_EVALUE))
        scored = []
        for group in loci(tert):
            score = sum(o[4] for _, o in group)
            cover = covered([iv for _, o in group for iv in o[5]]) / hmm_len
            has_rbd = any(name in rbd for name, _ in group)
            seq = group[0][1][0]
            locus = f"{seq}:{group[0][1][1]}-{max(o[2] for _, o in group)}:{group[0][1][3]}"
            scored.append((score, cover, len(group), has_rbd, locus))
        if scored:
            best = max(scored)
            row.update(best_score=f"{best[0]:.0f}", best_cover=f"{best[1]:.2f}",
                       best_orfs=best[2], best_rbd="yes" if best[3] else "no", best_locus=best[4])
            with_rbd = [x for x in scored if x[3]]
            if with_rbd:
                b = max(with_rbd)
                row.update(rbd_score=f"{b[0]:.0f}", rbd_cover=f"{b[1]:.2f}", rbd_locus=b[4])
        rows.append(row)
    with open(ROOT / "outputs" / "tert_hmm_summary.tsv", "w") as fh:
        print("\t".join(cols), file=fh)
        for r in rows:
            print("\t".join(str(r.get(c, "")) for c in cols), file=fh)
    print(f"[summarise_tert_hmm] {len(rows)} species -> outputs/tert_hmm_summary.tsv")


if __name__ == "__main__":
    main()
