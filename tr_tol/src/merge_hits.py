#!/usr/bin/env python3
"""Append each species' hits from earlier rounds to a round's hits.fa, so
the copy-aware template check sees every TR-like copy found so far.
Usage: merge_hits.py <round dir> <earlier round dir>..."""
import sys
from pathlib import Path
new, earlier = Path(sys.argv[1]), [Path(d) for d in sys.argv[2:]]
n = 0
for f in (new / "outputs").glob("*.hits.fa"):
    with open(f, "a") as out:
        for d in earlier:
            g = d / "outputs" / f.name
            if g.exists():
                out.write(g.read_text())
    n += 1
print(f"[merge_hits] {n} species")
