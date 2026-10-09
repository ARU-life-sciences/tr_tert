# De novo TR discovery by template co-variation

`src/covariation.py` finds TR candidates without a TR model, promoter or
structure: where the telomeric repeat differs between two related
species, the orthologous TR locus must carry a template matching each
species' own repeat. Chance template-like sequence has no reason to
change with the telomere.

1. Species A: every window of len(unit_A) + 2 nt exactly matching the
   tandem of unit_A (either strand), outside longer arrays (<= 3 copies),
   with 150 nt flanks.
2. blastn against species B (E <= 1e-10).
3. Co-varying: B's aligned window (+/- 8 nt) matches unit_B and not unit_A.

## Benchmarks against published Hymenoptera TRs (2026-10-09)

| A (unit) | B (unit) | template-like in A | aligned in B | co-varying | known TR |
|---|---|---|---|---|---|
| *Bombus sylvestris* (AACCCGAACCT) | *B. pratorum* (AACCCCAACCT) | 180 | 56 | **2** | **rank 1** |
| *Cerceris rybyensis* (AACCCAGACCT) | *Mimumesa dahlbomi* (AACCT) | 559 | 27 | 2 | not recovered |

Within a genus it works: two candidates genome-wide, the published TR
first, in about a minute. Between crabronid subfamilies it fails, and not
because of the filter: the two species' known TRs (from
`../tr_hym`) have no blastn similarity at all (best E = 0.24). TR sequence
outside the template turns over too fast for pairwise alignment beyond
close relatives, so the method needs a repeat change between congeners
(or near), as Fajkus et al. used *B. sylvestris* against other *Bombus*.

## Where it applies

27 genera in the ToL-wide run have >= 2 strongly called units, e.g.
*Polydrusus*, *Sitona*, *Otiorhynchus*, *Phyllobius* (weevils),
*Bruchidius*, *Agonum*, *Nomada*, *Crossocerus*, *Lasioglossum*, *Tiphia*,
*Biston*, *Catocala*, and plants (*Euphorbia*, *Lythrum*, *Ranunculus*,
*Hypericum*, *Lysimachia*, *Potentilla*). A TR found this way seeds a
profile model (HMM/CM) to search the rest of its clade, as the
Hymenoptera models do in `../tr_hym`.

## The 84-pair run (2026-10-09)

`src/make_pairs.py` -> `pairs.tsv` (84 pairs, 23 genera), `src/worker_pair.bash`
-> `outputs/pairs/<A>__<B>/{covarying.tsv,summary.txt}`; 73 finished.
- *Bombus sylvestris* against each of three congeners: 1 candidate each
  time, the published TR.
- Elsewhere tens to thousands of candidates, and the known TR (where
  `../tr_hym` has one) mostly not recovered: *Nomada* ranks 6 and 13,
  *Andrena* and *Lasioglossum* not at all. Requiring co-variation against
  every partner with another unit cuts the lists (*Nomada* 5-33) but lost
  the known TR in *N. goodeniana*. Pairs of short or one-base-different
  units (*Hypericum*, *Lythrum*, *Mechanitis*, *Tiphia*) are dominated by
  chance matches.
- The short-unit pairs wrote tens of GB each (218 GB in all) before the
  script deleted its intermediates; 11 unfinished pairs were stopped and
  removed, including six *Nomada flava*/*panzeri* pairs worth rerunning.

## Synteny-anchored search (2026-10-09)

`src/synteny_test.py`: 2 kb chunks from 100 kb either side of species A's
known TR mapped to B (minimap2 asm20); the anchors nearest A's TR give an
interval in B, searched for B's template. Tested on pairs with both TRs
known (`src/run_synteny_tests.bash` -> `outputs/synteny_tests*.tsv`):

| clade | B's TR inside the interval | interval |
|---|---|---|
| Hymenoptera | genus (*Bombus*) only; 0/7 beyond | 7 kb, 1 window = the TR |
| Lepidoptera | 7/9: genus to superfamily (butterflies vs *Noctua*); not Ditrysia-deep (too few anchors) | 8-10 kb (genus), 30-60 kb (family+) |

- In Lepidoptera the TR stays in place across families (stable
  chromosomes), so synteny cuts a genome to tens of kb. But the interval
  still holds 35-150 AACCT-like windows, so a second signal is needed
  inside it (local alignment of the TR body, or a model).
- In crabronids it fails for a biological reason: *Cerceris* and
  *Mimumesa* each have one validated TR, but each sits in a neighbourhood
  whose orthologue in the other species has none (both directions
  checked: flanks map cleanly to *Mimumesa* SUPER_15 / *Cerceris*
  SUPER_13). The TR gene has moved between these subfamilies.
