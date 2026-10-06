#!/usr/bin/env python3
"""Compare TR-paralog-differing sets with the TR-identical control group
on telomere array structure (notes/tr_core_template_array_structure.md,
headline #6), with two-sided Fisher's exact tests. Previously computed
by hand; scripted for re-running the control group on tidk 0.3.0 output
(notes/tidk_v03_rerun.md).

Reads two tr_core_template_array_structure_summary.py outputs (the
`structure` column):
  - regular alternation: REGULAR_ALTERNATION
  - any confirmed minor variant: anything but ABSENT
  - non-random structure, given present: BLOCKY, CLUSTERED or
    REGULAR_ALTERNATION, out of sets with a minor variant

Usage:
  array_structure_control_comparison.py <differing_summary.tsv> <control_summary.tsv>
"""
import csv
import sys

from scipy.stats import fisher_exact


def structures(path):
    return [r["structure"].split(" ")[0] for r in csv.DictReader(open(path), delimiter="\t")]


def compare(label, a, n_a, b, n_b):
    _, p = fisher_exact([[a, n_a - a], [b, n_b - b]])
    print(f"{label}\t{a}/{n_a} ({a / n_a:.1%})\t{b}/{n_b} ({b / n_b:.1%})\t{p:.3g}")


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    diff, ctrl = structures(sys.argv[1]), structures(sys.argv[2])
    nonrandom = {"BLOCKY", "CLUSTERED", "REGULAR_ALTERNATION"}
    print(f"comparison\tdiffering ({len(diff)})\tcontrol ({len(ctrl)})\tFisher p")
    compare("regular alternation", diff.count("REGULAR_ALTERNATION"), len(diff),
            ctrl.count("REGULAR_ALTERNATION"), len(ctrl))
    pd = [s for s in diff if s != "ABSENT"]
    pc = [s for s in ctrl if s != "ABSENT"]
    compare("any confirmed minor variant", len(pd), len(diff), len(pc), len(ctrl))
    compare("non-random, given present", sum(s in nonrandom for s in pd), len(pd),
            sum(s in nonrandom for s in pc), len(pc))


if __name__ == "__main__":
    main()
