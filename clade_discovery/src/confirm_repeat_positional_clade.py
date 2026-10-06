#!/usr/bin/env python3
"""Clade-discovery adaptation of the parent pipeline's
confirm_repeat_positional.py: tidk explore's top hit is the most
abundant tandem repeat GENOME-WIDE, which is not the same thing as
telomeric - satellite DNA, centromeric repeats, and other high-copy
tandem families routinely outrank the true (much shorter, array-bound)
telomeric repeat in a raw genome-wide count. This runs the actual
positional check (`tidk search`, restricted to chromosome termini) for
each species' top-N Stage-0 candidates, not just the single top hit,
since the true telomere motif isn't guaranteed to be #1 by raw count.

Usage:
  confirm_repeat_positional_clade.py <clade_name> <species> [top_n] [window_bp]

Reads candidates from clade_discovery/outputs/tidk/<species>.tidk.tsv
(Stage 0 output) and the genome path from
clade_discovery/inputs/<clade_name>_paths.txt.

Writes clade_discovery/outputs/positional/<species>_<repeat>_telomeric_repeat_windows.tsv
(raw tidk search output, one per candidate) and a per-species summary
to clade_discovery/outputs/positional/summaries/<species>.tsv - kept
one-file-per-species (rather than one shared appended file) so this can
be safely run as one LSF job per species in parallel with no risk of
concurrent-write corruption; merge_positional_confirmations.bash
concatenates them afterward.
"""
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

CLADE_ROOT = Path(__file__).resolve().parent.parent
TIDK = "tidk"
DEFAULT_TOP_N = 3
DEFAULT_WINDOW = 5000
MIN_TERMINAL_COUNT = 10
MIN_FOLD_OVER_MEDIAN = 5
SUMMARY_HEADER = ["clade", "species", "repeat", "explore_rank", "explore_count",
                   "n_main_chroms", "n_terminal_enriched", "verdict"]


def load_genome_path(clade_name, species):
    paths_file = CLADE_ROOT / "inputs" / f"{clade_name}_paths.txt"
    with open(paths_file) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == species:
                return parts[1]
    return None


def load_top_candidates(species, top_n):
    # tidk 0.3.0 output (<species>.explore.tsv, four columns) if present,
    # else tidk 0.2.7 (<species>.tidk.tsv, two columns) - see ../src/tidk_io.py
    tidk_dir = CLADE_ROOT / "outputs" / "tidk"
    tsv = tidk_dir / f"{species}.explore.tsv"
    if not tsv.exists():
        tsv = tidk_dir / f"{species}.tidk.tsv"
    if not tsv.exists() or tsv.stat().st_size == 0:
        return []
    rows = []
    with open(tsv) as fh:
        next(fh, None)  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1].isdigit():
                rows.append((parts[0], int(parts[1])))
    rows.sort(key=lambda x: -x[1])
    return rows[:top_n]


def run_tidk_search(species, repeat, genome_path, outdir, window):
    outdir.mkdir(parents=True, exist_ok=True)
    label = f"{species}_{repeat}"
    subprocess.run(
        [TIDK, "search", "-s", repeat, "-w", str(window), "-o", label,
         "-d", str(outdir), "--log", genome_path],
        check=True, capture_output=True,
    )
    return outdir / f"{label}_telomeric_repeat_windows.tsv"


def analyse(tsv_path):
    data = defaultdict(list)
    with open(tsv_path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            data[row["id"]].append((int(row["window"]),
                                     int(row["forward_repeat_number"]),
                                     int(row["reverse_repeat_number"])))

    MIN_WINDOWS = 20  # >=100kb at the default 5kb window - lower than the
                       # parent plant script's 100, since several of these
                       # clade genomes are much smaller than typical plant
                       # genomes (some beetle/spider assemblies are well
                       # under 500Mb) and would otherwise lose most/all
                       # scaffolds as "too short to count as a chromosome"
    MAX_MAIN_CHROMS = 40
    sized = sorted(data.items(), key=lambda kv: -len(kv[1]))
    main_ids = {chrom for chrom, windows in sized
                if len(windows) >= MIN_WINDOWS}
    if len(main_ids) > MAX_MAIN_CHROMS:
        main_ids = {chrom for chrom, _ in sized[:MAX_MAIN_CHROMS]}

    n_enriched = 0
    for chrom, windows in data.items():
        if chrom not in main_ids:
            continue
        windows = sorted(windows)
        totals = [f + r for _, f, r in windows]
        nonzero = [t for t in totals if t > 0]
        median = sorted(nonzero)[len(nonzero) // 2] if nonzero else 0
        first_total, last_total = totals[0], totals[-1]

        def enriched(count):
            return count >= MIN_TERMINAL_COUNT and (median == 0 or count >= MIN_FOLD_OVER_MEDIAN * median)

        if enriched(first_total) or enriched(last_total):
            n_enriched += 1

    return len(main_ids), n_enriched


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    clade_name = sys.argv[1]
    species = sys.argv[2]
    top_n = int(sys.argv[3]) if len(sys.argv) > 3 else DEFAULT_TOP_N
    window = int(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_WINDOW

    genome_path = load_genome_path(clade_name, species)
    if genome_path is None:
        sys.exit(f"species {species!r} not found in {clade_name}_paths.txt")

    candidates = load_top_candidates(species, top_n)
    if not candidates:
        print(f"[confirm_repeat_positional_clade] {species}: no Stage 0 candidates found, skipping", file=sys.stderr)
        return

    outdir = CLADE_ROOT / "outputs" / "positional"
    summaries_dir = CLADE_ROOT / "outputs" / "positional" / "summaries"
    summaries_dir.mkdir(parents=True, exist_ok=True)
    summary_tsv = summaries_dir / f"{species}.tsv"
    with open(summary_tsv, "w") as out:
        out.write("\t".join(SUMMARY_HEADER) + "\n")
        for rank, (repeat, count) in enumerate(candidates, 1):
            print(f"[confirm_repeat_positional_clade] {species}: rank {rank} "
                  f"'{repeat}' (explore count {count})", file=sys.stderr)
            tsv_path = None
            try:
                tsv_path = run_tidk_search(species, repeat, genome_path, outdir, window)
                n_chroms, n_enriched = analyse(tsv_path)
            except subprocess.CalledProcessError as e:
                print(f"  ERROR running tidk search: {e}", file=sys.stderr)
                continue
            finally:
                # Raw per-window tidk search output is 10-20MB+ per
                # species/candidate - at ~190 species x top_n candidates
                # this adds up fast, and this farm's disk quota has
                # already been hit once this project (see
                # notes/organelle_and_cobiont_assembly_pitfall.md's
                # sibling incident in the parent repo's session history).
                # The summary row captures everything this pipeline
                # actually needs; delete the raw dump rather than
                # accumulate it - regenerate with run_tidk_search() above
                # if a specific case needs deeper inspection.
                if tsv_path is not None:
                    Path(tsv_path).unlink(missing_ok=True)
            verdict = ("TERMINAL_ENRICHED" if n_enriched >= 2 else
                       "WEAK" if n_enriched == 1 else "NOT_TERMINAL")
            out.write("\t".join(str(x) for x in [
                clade_name, species, repeat, rank, count, n_chroms, n_enriched, verdict
            ]) + "\n")
            print(f"  -> {n_enriched}/{n_chroms} chromosomes terminal-enriched: {verdict}", file=sys.stderr)


if __name__ == "__main__":
    main()
