# Paper outline: telomere evolution across the Tree of Life (anchor, 2026-10-10)

The plan for the phase 4 paper (tidk roadmap). Background and prior work:
`literature_synthesis.md`. Working title: *How telomeres evolve: repeat
changes, telomerase RNA mobility and telomerase loss across 3,400 genomes.*

## Thesis

Telomeric repeats change by two routes - single template mutations, and
larger changes that keep the ancestral repeat and add to it - under a
constraint that every repeat keeps a G-run. The large changes are made
possible by telomerase RNA (TR) duplication and relocation, which is frequent
in some lineages (Hymenoptera) and rare in others (Lepidoptera). Where
telomerase is lost altogether, other things occupy the chromosome ends.

## 1. Repeat evolution (the backbone)

Data: tidk 0.3.0 calls checked at chromosome ends, 3,422 DToL species
(`a-telomeric-repeat-database`, `tol_run/`).

- **Two routes of change.** Of 237 strongly called non-ancestral repeats (by
  species), 96 are one edit from their group's ancestral repeat (`AACCT` →
  `ACCTG`, `AACCC`, `AAACC`; plant `AAACCCT` → `AACCCT`), 123 are four or more
  (`AACCCAGACCT`, `AACCATCCCT`, `AAAGAACCT`, ...), and only 18 lie between.
  Most large changes contain the ancestral repeat plus an insertion.
  To do: per-lineage (not per-species) counts; test against TR copy state.
- **A G-run constraint.** All 268 non-ancestral repeats keep `GG` or `GGG` on
  the G-rich strand (187 `GGG`). To do: compare with random units of the same
  length and composition; G-quadruplex potential; telomere-binding protein
  (TRF/TRB, POT1) coevolution.
- **Rates and patterns on the tree.** Number of independent changes,
  reversals, direction, repeated targets (`AACCC`, `ACCTG` in several beetle
  and wasp lineages). Needs a ToL species tree.
- **Mixed and compound telomeres** as lineages in transition (two working
  templates).
- **Predictors:** holocentric chromosomes (Hemiptera diverse vs Lepidoptera
  uniform), genome size, GC, karyotype change, TERT/TR copy state.

## 2. Mechanism: telomerase RNA

Data: `tr_tert/tol/` (TR.md, tr_hym, denovo, tr_copies, tr_tol).

- **Templates validate the repeats:** where a TR is found, its template
  encodes the tidk repeat in 1,087/1,096 strongly called species.
- **TRs move:** positional survey, 250 pairs - conserved in Lepidoptera (genus
  25/25, family 22/25), not in Hymenoptera (16/25, 4/24); orthologous interval
  conserved and empty. Precedents in one genus each: *Mimulus* (Kumawat 2023),
  *Andrena* (Klapproth 2025). To do: extend to vertebrates, plants, new animal
  groups; predictors of mobility.
- **TRs duplicate:** Hymenoptera mean 3.3 TR-like copies vs Lepidoptera 1.2
  (p = 3e-22); pimplines 3-7 copies, one with a template; halictids up to 52.
- **Two routes explained?** Test: large repeat changes coincide with
  multi-copy or relocated TRs, small ones with single-copy TRs.
- **Mechanism of movement:** no retrocopy signature; dispersed TR-sized
  insertions vs copies inside larger repeats. To do: transposon annotation
  (Kumawat's DNA-transposon hypothesis), copy ages.
- **TERT coevolution** (*Mimulus* yes, *Andrena* no): rates of TERT change on
  branches with TR duplication, relocation or a repeat shift; ~300 TERT
  proteins in hand.
- **Biogenesis switches:** classify TRs (type-3 promoter, box H/ACA,
  intronic) and map them.

## 3. Telomerase loss

- TERT absent and no telomere at chromosome ends in Diptera (expected) and
  spiders (0/6 TERT, 0/18 telomeres: new), *Viscum*.
- To do: map losses tree-wide; characterise chromosome ends without
  telomerase (retrotransposon arrays, satellites, recombination signatures);
  interrupted telomeres (Lepidoptera `not_terminal`, TRAS/SART) as a possible
  first step.

## 4. Further results (supporting or separate)

- Telomere length atlas: array lengths at T2T ends (`tidk ends`) and per-read
  lengths from HiFi (`tidk length`), against body size, lifespan, life history.
- Interstitial telomeric sequences as dated fusion scars (Lepidoptera Merian
  elements).
- TR horizontal transfer: the weevil TR 96% identical to a riodinid butterfly
  TR; viral precedent (Marek's disease virus).
- T2T status of every DToL assembly as a QC resource.

## Methods (from phases 1-3 of the tidk roadmap)

tidk 0.3.0 (`explore` with primitive units and error tolerance, `ends` for
chromosome-end checks, `length`); TR search by model walking with
telomere-validated seeds, lineage-specific models, copy-aware template
checks, covariance models, synteny and co-variation; TERT by miniprot + RBD
domain and profile HMMs.

## Status

Done: 1 (data, first numbers), 2 (validation, relocation survey, copies), 3
(TERT survey). Running: ToL-wide TR search with literature seeds (molluscs,
cnidarians, echinoderms, plants beyond flowering plants, fungi). Missing: a
species tree; per-lineage counts; the tests marked "To do".
