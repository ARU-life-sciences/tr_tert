#!/usr/bin/env python3
"""ToL-wide TR model walking, one round.

Seeds: every validated TR from earlier rounds (template_check.tsv, run
copy-aware: strong telomere call, template supports it, E <= 1e-10). The
seed sequence is the hit at the template_check locus - not simply the
lowest-E hit, which in multi-copy genomes is often a decayed copy.

Models: seeds grouped lineage by lineage - by family where a family has
>= MIN_GROUP seeds, otherwise by order - so models stay specific (a pooled
Ichneumonidae model missed the pimplines; sawfly TRs filed under the wrong
group were invisible). The round searches only these seed models; the
stage-1 library has already run on every genome, and its hits are merged
in before the template check (src/merge_hits.py).

Targets: every species without a validated TR.

Usage: walk_tol.py <new round dir> <earlier round dir>...
  e.g. walk_tol.py ../tr_tol_r2 ../tr_tol
Then: CLADE_DIR=<new> MODEL_NAME=seeds bash ../tr_hym/src/build_models.bash,
cat into models/library.hmm (done here if build succeeds), and run
src/worker_tol.bash with CLADE_DIR=<new>.
"""
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

MAX_EVALUE = 1e-10
MIN_GROUP = 3
FLANK = 300
HERE = Path(__file__).resolve().parent.parent  # tr_tol


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].strip(), []
        else:
            seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def core(header, seq):
    loc = header.split("|")[0]
    a, b = (int(x) for x in loc.rsplit(":", 1)[1].split("-"))
    off = a - max(1, a - FLANK) if a < b else FLANK
    return seq[off:off + abs(b - a) + 1].upper()


def validated(round_dir):
    """species -> TR sequence at the template_check locus."""
    out = {}
    for r in csv.DictReader(open(round_dir / "outputs" / "template_check.tsv"), delimiter="\t"):
        if (r["agrees"] != "yes" or r["tidk_support"] != "strong" or not r["tr_evalue"]
                or float(r["tr_evalue"]) > MAX_EVALUE):
            continue
        for h, s in read_fasta(round_dir / "outputs" / f"{r['species']}.hits.fa"):
            if h.split("|")[0] == r["tr_locus"]:
                out[r["species"]] = core(h, s)
                break
    return out


def main():
    new = Path(sys.argv[1]).resolve()
    earlier = [Path(d).resolve() for d in sys.argv[2:]]
    species = list(csv.DictReader(open(HERE / "species.tsv"), delimiter="\t"))
    info = {r["species"]: r for r in species}
    seeds = {}
    for d in earlier:
        seeds.update(validated(d))
    fam = defaultdict(list)
    for sp, s in seeds.items():
        fam[(info[sp]["order"] or info[sp]["class"] or "NA", info[sp]["family"] or "NA")].append((sp, s))
    groups = defaultdict(list)
    for (order, family), items in fam.items():
        key = f"{order}_{family}" if len(items) >= MIN_GROUP else order
        groups[key] += items
    for d in ("models", "outputs", "logs"):
        (new / d).mkdir(parents=True, exist_ok=True)
    for g, items in groups.items():
        with open(new / "models" / f"{g}.fa", "w") as fh:
            for sp, s in items:
                fh.write(f">{sp}\n{s}\n")
    todo = [r for r in species if r["species"] not in seeds]
    with open(new / "species.tsv", "w") as fh:
        print("\t".join(species[0].keys()), file=fh)
        for r in todo:
            print("\t".join(r.values()), file=fh)
    print(f"[walk_tol] {len(seeds)} seeds in {len(groups)} groups; {len(todo)} species to search -> {new}")
    # build the seed HMMs and the round's library: stage-1 library + seeds
    subprocess.run(["bash", str(HERE.parent / "tr_hym" / "src" / "build_models.bash")],
                   env={"CLADE_DIR": str(new), "MODEL_NAME": "seeds", "PATH": "/usr/bin:/bin"},
                   check=True, stdout=subprocess.DEVNULL)
    # the round searches only the new seed models: the stage-1 library (HMMs
    # and CMs) has already been run on every genome, and its hits are
    # merged back in before the template check (src/merge_hits.py)
    with open(new / "models" / "library.hmm", "w") as out:
        for line in open(new / "models" / "seeds.hmm"):
            out.write(line.replace("NAME  TR_", "NAME  seed_TR_"))


if __name__ == "__main__":
    main()
