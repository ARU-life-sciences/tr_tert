# TR in Hymenoptera: does the template encode the telomere? (2026-10-09)

Published Hymenoptera TRs (Fajkus et al. 2023, NAR, 10.1093/nar/gkac1202:
256 TRs, Table S1 GenBank accessions; supplement in
`../inputs/fajkus2023/`, fetched from Europe PMC) -> per-family DNA
profile HMMs (`src/build_models.bash`, 14 models) -> nhmmer against all
298 Hymenoptera in the ToL-wide run (`src/worker_search.bash`) -> the
template in each species' best hit (E <= 1e-10) tested against the tidk
call (`src/template_check.py` -> `outputs/template_check.tsv`).

A template supports a unit if it exactly matches >= len(unit) + 2 nt of
its tandem (templates are often ~1.2 units: *Bombus*, 13 nt of
AACCCCAACCT).

| tidk call | TR found, template = tidk call | TR found, other | no TR hit |
|---|---|---|---|
| strong (>= 25% of ends) | **110** | 5 | 74 |
| weak | 26 | 11 | 31 |

- **Where a TR is found, its template encodes the repeat at the chromosome
  ends in 110/115 strongly called species**, across every shift lineage:
  AACCCAGACCT (24), AAACAAACCC (17), AAAGAACCT, AACCCCAACCT, AACCCAGACGC,
  AAACCTAACCC, ... The 5 others: *Agrypon varitarsum* (tidk called the
  doubled unit, AAACAAACCCAAACCC; the template is AAACAAACCC), three
  *Andrena* with unusual calls (the genus has 7 different units in 8
  species), and *Alloplasta piceator* (no template-like match).
- The published repeats agree with the tidk calls where they overlap:
  TTAGGTCTGGG = AACCCAGACCT, TTAGGTTGGGG = AACCCCAACCT, TTTGTTTGGG =
  AAACAAACCC.
- No TR hit: mostly families with few or no published TRs - sawflies
  (Tenthredinidae 7, Cephidae 4, Argidae 3, Athaliidae 2, ...),
  Ichneumonidae (31 of ~60; 22 published TRs), Vespidae, Chrysididae.
  Next: add the found DToL TRs to the models and search again, as for TERT.
- *Andrena denticulata* hung in nhmmer twice and was stopped (not run).
- nhmmer needs `--dna`: T2T assemblies can start with kilobases of G-free
  telomere, which defeats its alphabet guess (79 species failed without it).
