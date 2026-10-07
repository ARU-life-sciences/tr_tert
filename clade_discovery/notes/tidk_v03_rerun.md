# Clade discovery: tidk 0.3.0 rerun of Stage 0 / 0.5

**Status: done for all 188 species (2026-10-07).** Results in
`outputs/tidk_v0.3/` (explore + `summary.tsv`) and `outputs/tidk_ends/`.
The v0.2.7 Stage 0 results in `outputs/tidk/` and the Stage 0.5
`positional_confirmations.tsv` are left untouched.

Same method as the parent project's rerun (`../notes/tidk_v03_rerun.md`):
tidk 0.3.0 `explore` with the Stage 0 parameters, then `tidk ends` for each
species' top 5 explore units, calling the telomere as the candidate found
at the most chromosome ends (sequences >= 1 Mb, 10 kb windows). This
replaces Stage 0 + 0.5 in one pass, checking 5 candidates rather than 3,
from rankings without tidk 0.2.7's multimer and counting bugs. Run with
the parent `src/run_tidk_v03.bash`, `run_tidk_ends_candidates.bash` and
`tidk_v03_summary.py`, pointed at `inputs/{coleoptera,araneae,myriapoda}_paths.txt`.

## Results

| clade | species | telomere at >= 1 chromosome end | T2T chromosomes | Stage 0.5 (v0.2.7) |
|---|---|---|---|---|
| Coleoptera | 166 | **159** | 515/3017 | 106 |
| Araneae | 18 | **0** (best candidate at most 1 end) | 0/277 | 4 |
| Myriapoda | 4 | 4 (`AACCT`) | 3/51 | 1 |

In 51 beetles the telomere is not the top explore unit (rank 2-5), the
same satellite-outranks-telomere effect seen in plants.

**Coleoptera telomere units**: `AACCT` (`TTAGG`, the ancestral insect
repeat) in 104 species, plus three variant families, all with balanced
strands (a real telomere reads both ways round):

- **`AACCC`** in *Aphodius fimetarius*, *A. granarius*, *A. sticticus*,
  *Geotrupes spiniger* (scarabs) and *Oxystoma pomonae* (a weevil): 6-16
  chromosome ends each, minor strand 28-49%. These are the Stage 0.5
  "minority variant" species, plus *A. fimetarius*. *A. ater* (a Stage 0.5
  `AACCC` species) is ambiguous: `AACCC` is still present (1196 copies,
  balanced) and at 6/62 ends, but it is explore rank 24, outside the top 5
  candidates; the call is the top unit `AAACGATAAT`, at 15/62 ends, which
  may be a satellite immediately next to the telomere rather than the
  telomere itself. All five assemblies are PacBio HiFi + Hi-C with no ONT
  data, so this is not the ONT basecalling artefact (`TTAGG` read as
  `TTGGG`, i.e. `AACCC`) seen in 2019 ONT reads; together with its
  clustering in scarabs, it looks like a real telomere variant.
- **`ACCTG`** (`TCAGG`, the *Tribolium* repeat) in 16 species, well
  supported in e.g. *Anthocomus fasciatus* (15/20 ends), *Cis
  bilamellatus*, *Crepidodera* spp., *Mycetophagus* spp. (8-11 ends), weak
  in several *Anaspis* (1-7/16) and *Quedius lateralis* (0/36).
- **`ACCCT`** in *Carabus intricatus*, *C. nemoralis*, *C. problematicus*
  (18-21 of ~30 ends each).

**Araneae: no telomere at any chromosome end in any of the 18 species.**
4 have no tandem repeat near their ends at all; every other call is a
satellite or low-complexity repeat at 0-1 ends. This resolves the Stage
0.5 open item (its 4 "confirmed" spider motifs looked like low-complexity
artefacts): none of them is at the chromosome ends. The same pipeline
finds telomeres in 159/166 beetles, so this is not a tool failure. It is
consistent with spiders lacking the arthropod `TTAGG` telomere, and with
Fajkus & Fajkus 2025 listing spiders as a telomerase-loss candidate, but
assembly-end completeness is the caveat: chromosome ends that are
missing from the assembly would look the same.

**Caveats of the telomere call** (both seen in *A. ater*): only the top 5
explore units are checked with `tidk ends`, so a telomere ranked lower is
missed; and "found at the most chromosome ends" can pick a satellite that
sits immediately next to the telomere (within the 1 kb `--max-offset`).
Calls at only a few ends deserve a manual look.

## Implications for the discovery pipeline

- The `AACCC` scarab group is the natural seed for reciprocal filtering
  (Method step 4), as already planned - now confirmed at the chromosome
  ends with strand balance, and with 4 species rather than relying on
  top-3 positional checks.
- `ACCTG`-family beetles are a second, larger variant group (16 species
  across several families) that wasn't visible in the v0.2.7 rankings.
- For Araneae, telomere-motif-based TR discovery has no motif to work
  with; if spiders have lost telomerase, the TERT/TR-absence angle
  (as for `Viscum` in the parent project) is the relevant one.
