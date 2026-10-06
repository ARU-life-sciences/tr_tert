"""Reading `tidk explore` output from either tidk version used in this project.

v0.2.7 - outputs/tidk/<species>.tidk.tsv (src/run_tidk.bash):
    canonical_repeat_unit  count_repeat_runs_gt_<threshold>
  Repeats are the tandemly repeated kmer seeds themselves, so the list
  includes multimers and compound strings (e.g. AAACCCTAAACCCTAAACCCT).

0.3.0 - outputs/tidk_v0.3/<species>.explore.tsv (src/run_tidk_v03.bash):
    canonical_repeat_unit  copies  copies_as_unit  copies_as_revcomp
  Repeats are primitive units: exact multimers are collapsed into the unit,
  counts are copies of the unit, and runs tolerate sequencing errors.

The directory comes from the TIDK_DIR environment variable, defaulting to
outputs/tidk so existing results reproduce unchanged. Set
TIDK_DIR=outputs/tidk_v0.3 (absolute, or relative to the project root) to
use 0.3.0 results.
"""
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def tidk_dir():
    d = Path(os.environ.get("TIDK_DIR", "outputs/tidk"))
    return d if d.is_absolute() else PROJECT_ROOT / d


def tidk_path(species, directory=None):
    """The explore output for a species, preferring the 0.3.0 file name."""
    directory = Path(directory) if directory else tidk_dir()
    for name in (f"{species}.explore.tsv", f"{species}.tidk.tsv"):
        path = directory / name
        if path.exists():
            return path
    return directory / f"{species}.tidk.tsv"


def load_explore(path):
    """(version, [(repeat, count), ...]) in the order tidk ranked them.
    version is "0.2" or "0.3", from the header; ([], None) if missing."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return None, []
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        version = "0.3" if len(header) >= 2 and header[1] == "copies" else "0.2"
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1].isdigit():
                rows.append((parts[0], int(parts[1])))
    return version, rows


def tandemise(unit, length):
    """`length` bp of the unit repeated in tandem, e.g. to compare a
    primitive 0.3.0 unit with sequences longer than one copy."""
    return (unit * (length // len(unit) + 1))[:length]
