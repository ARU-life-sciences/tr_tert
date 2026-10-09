#!/usr/bin/env python3
"""A covariance model (Infernal) from a set of TR sequences.

HMM model walking (../tr_hym_r*) stops where sequence similarity runs out;
a CM also scores conserved secondary structure, which TRs keep long after
their sequence has diverged. Steps:
  1. MAFFT L-INS-i alignment of the sequences
  2. consensus secondary structure with RNAalifold (ViennaRNA), which
     scores covariation (compensatory base changes) across the alignment
  3. Stockholm alignment with #=GC SS_cons, then cmbuild
cmcalibrate (slow) is run separately, before cmsearch (src/run_cm.bash).

Usage: build_cm.py <seqs.fa> <name> <outdir>
Writes <outdir>/<name>.{aln.fa,sto,cm} and prints the consensus structure.
"""
import subprocess
import sys
from pathlib import Path

import RNA

MAFFT = "/software/team301/mafft-7.525-with-extensions/scripts/mafft"
MAFFT_ENV = {"MAFFT_BINARIES": "/software/team301/mafft-7.525-with-extensions/binaries",
             "PATH": "/usr/bin:/bin"}
CMBUILD = "/software/team301/infernal-1.1.4-linux-intel-gcc/binaries/cmbuild"


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def main():
    fa, name, outdir = sys.argv[1:4]
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    aln = subprocess.run([MAFFT, "--localpair", "--maxiterate", "1000", "--thread", "8",
                          "--quiet", fa], capture_output=True, text=True, check=True,
                         env=MAFFT_ENV).stdout
    (out / f"{name}.aln.fa").write_text(aln)
    tmp = out / f"{name}.aln.fa"
    seqs = [(n.replace("|", "_")[:60], s.upper().replace("T", "U")) for n, s in read_fasta(tmp)]
    structure, energy = RNA.alifold([s for _, s in seqs])
    width = max(len(n) for n, _ in seqs) + 2
    with open(out / f"{name}.sto", "w") as fh:
        fh.write("# STOCKHOLM 1.0\n\n")
        for n, s in seqs:
            fh.write(f"{n:<{width}}{s.replace('-', '.')}\n")
        fh.write(f"{'#=GC SS_cons':<{width}}{structure}\n//\n")
    subprocess.run([CMBUILD, "-F", "-n", name, out / f"{name}.cm", out / f"{name}.sto"],
                   check=True, stdout=subprocess.DEVNULL)
    pairs = structure.count("(")
    print(f"[build_cm] {name}: {len(seqs)} sequences, {len(seqs[0][1])} columns, "
          f"{pairs} consensus base pairs, alifold energy {energy:.1f}")
    print(structure)


if __name__ == "__main__":
    main()
