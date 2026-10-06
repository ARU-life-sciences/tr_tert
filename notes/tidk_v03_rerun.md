# tidk 0.3.0 rerun: corrected telomeric repeat calls for all plant assemblies

**Status: rerun complete for all 703 species (2026-10-07).** Results in `outputs/tidk_v0.3/` and `outputs/tidk_ends/`;
per-species summary in `outputs/tidk_v0.3/summary.tsv`. The v0.2.7
results in `outputs/tidk/` are left untouched, as the existing notes cite
them.

## Why rerun

Every `tidk explore` result in this project so far came from tidk 0.2.7,
which had bugs found and fixed in tidk 0.3.0 (see the tidk CHANGELOG):

- **Multimers reported as separate repeats.** `AAACCCTAAACCCTAAACCCT` and
  `AAACCCT` were separate rows, so the top "repeat" was often a multimer
  or a compound string rather than the unit.
- **The first run on each chromosome arm started at position 0**, so it
  counted everything before it. This inflated counts and let low-complexity
  junk (e.g. `AAAAAAAAAAG`) reach the top of the list.
- **Run counts double-counted some runs and dropped others**, and counted
  the same telomere again at every kmer length in a `-m/-x` range.

0.3.0 collapses repeats to their primitive canonical unit, merges runs
found at different kmer lengths, tolerates sequencing errors within runs,
and reports how many copies read as the unit vs its reverse complement
(`copies_as_unit` / `copies_as_revcomp`). It also adds `tidk ends`,
which calls a telomere (or not) at both ends of each chromosome.

## Method

- `src/run_tidk_v03.bash` / `worker_tidk_v03.bash`: `tidk explore` with
  the same parameters as the v0.2.7 run (`-m 4 -x 30 --threshold 100
  --distance 0.01`), then `tidk ends` with the top unit, one LSF job per
  species. Binary pinned at `bin/tidk-0.3.0` (tidk v0.3.0 tag).
- `src/run_tidk_ends_candidates.bash`: `tidk ends` for each species' top
  5 explore units. **The top explore unit isn't always the telomere**: a
  long subterminal satellite can have more copies near the chromosome
  ends than the telomere itself. So the telomere is called as the
  candidate found at the most chromosome ends (`src/tidk_v03_summary.py`).
- `tidk ends` only checks sequences of at least 1 Mb (no unplaced
  scaffolds or organelles), in 10 kb windows at each end.

## Results

- **The canonical land plant repeat (`AAACCCT`, i.e. `TTTAGGG`) is the
  telomere in 670/700 species.** Using v0.2.7's top repeat, even after
  collapsing it to its primitive unit, it was top in only 173/702.
  v0.2.7 and 0.3.0 top repeats compare as: 168 the same unit, 271 a
  compound/multimer string of the same unit, 261 different (many of these
  low-complexity junk at the top of the v0.2.7 list).
- **The top explore unit was not the telomere in 30 species** (telomere at
  explore rank 2-5), e.g. `Fragaria_vesca`, where a 12 bp satellite
  (`AGCCCCCTGCCG`) has more copies than `AAACCCT` but is at 0/14
  chromosome ends, while `AAACCCT` is at 10/14.
- **30 species have a non-canonical telomere call:**

| Species | Telomere unit | Explore rank | Ends with telomere |
|---|---|---|---|
| *Cymbalaria muralis* | `AAAAAACCCT` | 1 | 11/14 |
| *Anagallis arvensis* | `AAAACCCT` | 1 | 19/40 |
| *Saxifraga aizoides* | `AAAACCCT` | 1 | 18/26 |
| *Saxifraga cespitosa* | `AAAACCCT` | 1 | 35/86 |
| *Glaucium flavum* | `AAAACCCTACCCG` | 1 | 6/16 |
| *Chrysosplenium oppositifolium* | `AAAACCG` | 1 | 37/42 |
| *Hypericum androsaemum* | `AAAAGGCAAG` | 1 | 15/40 |
| *Hypericum hircinum* | `AAAAGGCAAGAAAAGGCATG` | 3 | 18/56 |
| *Mimulus peregrinus* | `AAACCCG` | 1 | 83/106 |
| *Ceratocapnos claviculata* | `AAACCG` | 1 | 24/32 |
| *Fumaria muralis* | `AAACCG` | 1 | 8/16 |
| *Carlina vulgaris* | `AAACCGAACCG` | 1 | 13/26 |
| *Potentilla anglica* | `AAACTCT` | 1 | 35/58 |
| *Crocus tommasinianus* | `AACCCT` | 1 | 8/24 |
| *Diapensia lapponica* | `AACCCT` | 1 | 10/14 |
| *Euphorbia amygdaloides* | `AACCCT` | 1 | 4/28 |
| *Euphorbia maculata* | `AACCCT` | 1 | 25/48 |
| *Euphorbia peplus* | `AACCCT` | 1 | 12/16 |
| *Helianthemum oelandicum* | `AACCCT` | 1 | 9/22 |
| *Iris foetidissima* | `AACCCT` | 1 | 32/116 |
| *Lythrum portula* | `AACCCT` | 1 | 8/10 |
| *Hypericum hirsutum* | `AACCCTG` | 1 | 5/26 |
| *Hypericum pulchrum* | `AACCCTG` | 1 | 12/18 |
| *Jasione montana* | `AACCCTG` | 1 | 0/16 |
| *Ranunculus repens* | `AACCCTG` | 2 | 21/36 |
| *Ranunculus sardous* | `AACCCTG` | 1 | 4/64 |
| *Solanum dulcamara* | `AACCCTG` | 1 | 11/24 |
| *Solanum nigrum* | `AACCCTG` | 1 | 26/72 |
| *Eleocharis palustris* | `AACTTGTACG` | 1 | 1/108 |
| *Calendula officinalis* | `ACTGCT` | 1 | 1/36 |

  Several of these look biologically coherent (e.g. `AACCCT` /
  `TTAGGG` in Asparagales and *Euphorbia*; `AAACCG`-type variants in
  Fumarioideae and *Chrysosplenium*; `AAACCCG` in *Mimulus peregrinus*,
  the genus Kumawat et al. 2025 studied for TR divergence; `AACCCTG` in
  *Solanum*, *Hypericum* and *Ranunculus*). Calls with almost no ends
  supporting them (*Eleocharis palustris* 1/108, *Calendula officinalis*
  1/36, *Jasione montana* 0/16, *Ranunculus sardous* 4/64) are not real
  telomere calls: no candidate is at the chromosome ends.
- **No tandem repeat near the chromosome ends at all** (explore found
  nothing above the threshold): `Hottonia_palustris`,
  `Sarmentypnum_sarmentosum`. Worth checking.
- **Telomere-to-telomere status:** 5584/12856 chromosome-sized sequences
  (43.4%) have a telomere at both ends; 12 assemblies are fully T2T. 41
  ends in 25 species have a telomere on the wrong strand (C-rich expected
  at sequence starts, G-rich at ends), suggesting misjoins or inverted
  ends, and 434 ends in 243 species have a telomere further than 1 kb from
  the end (`not_terminal`). This is per-assembly context for any claim
  that something is *absent* from a genome.

- **No telomere at any chromosome end** (no candidate at any end):
  `Impatiens_glandulifera`, `Jasione_montana`, `Ribes_uva_crispa`,
  `Senecio_squalidus`, `Tragopogon_pratensis` - and effectively
  `Viscum_album` (below). Probably assembly-related in most cases; worth a
  look before any of these is used in an absence argument.

## Impact on existing results

- **#1 Viscum:** no telomere at the chromosome ends. Canonical `AAACCCT` is
  at 0/676 sequence ends (338 sequences >= 1 Mb, 49 of them >= 100 Mb, the
  largest ~2.1 Gb), and the best of the top 5 explore candidates (a
  satellite, `AATAATTAGT`) is at 1/676. `AAACCCT` is present in the
  genome but sparse (explore rank 27, 1009 copies near ends). This is
  consistent with the missing TR and TERT, but only as supporting
  evidence: ~1% of assemblies here have no telomere at any end for
  probably assembly-related reasons, and the assembly-completeness gap
  noted for #1 still applies.
- **#2 Carlina:** canonical `AAACCCT` is at 0/26 chromosome ends, so the
  non-canonical claim stands. But the telomere is a mixture, not
  `AAACTG`-family alone: `AAACCGAACCG` (and `AACCG`) is at 13/26 ends
  and is the most abundant unit (27,083 copies); `AAACTG` is at 11/26
  ends (877 copies as a unit). The Template-to-repeat matches are
  unchanged (12 bp and 11 bp exact matches, shuffle p ~ 0, now against
  `AAACTG` and `AAACCGAACTGAAACTG`). Suggest describing it as a mixed
  `AAACCG`/`AAACTG`-family telomere.
- **#3/#4 and the template correlation:** `tr_repeat_correlation.py --all`
  on 0.3.0 output agrees with the v0.2.7 run at 1327/1430 loci (93%);
  significant matches (p < 0.05, >= 8 bp) 1070 -> 1141 (16 lost, spread
  thinly; 87 gained, 44 in `Levisticum_officinale`). Positional
  confirmations use `tidk search`, which these fixes don't affect.
- **#6 control group:** on 0.3.0 output,
  `find_identical_set_second_candidates.py` selects 60 control sets (62 on
  v0.2.7); 36 are shared, 26 of those with the same second candidate (up
  to rotation). (A first adapter tandemised every shorter unit to the core
  length, which gave 196 sets; that was wrong - a 0.2.7 seed of length L
  means an array with period L, so the 0.3.0 equivalent is only units whose
  length divides L, i.e. the multimers 0.2.7 reported.) The thresholds
  were left as they were: 0.3.0 counts are not a fixed multiple of 0.2.7's
  (0.3.0/0.2.7 for the `AAACCCT` family: median 0.39, IQR 0.09-2.97 across
  671 species, as 0.2.7's were inflated by the bugs), so there is no
  principled rescaling. Re-running the base-pair terminal walk on the
  0.3.0 set (`src/run_identical_set_array_structure_bsub.bash`, one LSF
  job per species; outputs in `outputs/tr_repeat_correlation/tidk_v03/`;
  Fisher tests scripted in `src/array_structure_control_comparison.py`,
  which reproduces the published v0.2.7 numbers exactly):

  | | Differing (44) | Control, v0.2.7 (62) | Control, 0.3.0 (60) |
  |---|---|---|---|
  | Regular alternation (candidate HOR) | 6/44 (13.6%) | 0/62, p=0.004 | **0/60, p=0.005** |
  | Any confirmed minor variant | 41/44 (93.2%) | 45/62 (72.6%), p=0.011 | **55/60 (91.7%), p=1.0** |
  | Non-random structure, given present | 35/41 (85.4%) | 30/45 (66.7%), p=0.049 | 44/55 (80.0%), p=0.59 |

  At `GAP_THRESHOLD=100`: regular alternation 6/44 vs 0/60 (p=0.005); any
  minor variant 95.5% vs 95.0% (p=1.0); non-random given present 85.7% vs
  87.7% (p=0.77).

  **The headline holds**: regular alternation is still found only in
  TR-paralog-differing sets, at both gap thresholds. **The presence
  result does not**: with 0.3.0 candidates, a minor variant is confirmed
  in as many controls as differing sets. Most likely the v0.2.7
  candidates were partly bug artefacts not genuinely in the telomere,
  which made controls look as though they lacked minor variants - so
  that comparison should be dropped. The third comparison was already
  not robust and stays non-significant.

## Scripts changed to read 0.3.0 output

`src/tidk_io.py` reads either version (detected from the header) and
finds either file name; `TIDK_DIR` selects the directory, defaulting to
`outputs/tidk` so existing results reproduce unchanged (checked
byte-identical for `tr_repeat_correlation.py` and
`find_identical_set_second_candidates.py`). Set
`TIDK_DIR=outputs/tidk_v0.3` for 0.3.0. Where matching relied on 0.2.7's
tandemised strings, 0.3.0 units are tandemised back
(`tr_repeat_correlation.py`, `find_identical_set_second_candidates.py`).
`clade_discovery/src/confirm_repeat_positional_clade.py` also accepts both.
Before these changes, `tr_repeat_correlation.py` and
`confirm_repeat_positional_clade.py` would have read zero repeats from
0.3.0 output, silently.

## Next

- Rerun the `clade_discovery` species (188) the same way.
- A Tree of Life-wide run (all taxa) with the same scripts and an
  all-taxa assembly list.
