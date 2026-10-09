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
