# Paper outline: telomere evolution across the Tree of Life (anchor, 2026-10-10)

The plan for the phase 4 paper (tidk roadmap). Background and prior work:
`literature_synthesis.md`. Working title: *How telomeres evolve: repeat
changes, telomerase RNA mobility and telomerase loss across 3,400 genomes.*

## Thesis

Telomeric repeats change by two routes - single template mutations, and
larger changes that keep the ancestral repeat and add to it - with
G-triplets strongly enriched (a G-quadruplex or protein-binding constraint). The large changes are made
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
  Caveat: these counts only included units containing `CC` (as the first
  G-run count did). **Per lineage** (B1/B2, `repeat_evolution/transitions.py`:
  Fitch parsimony on the Open Tree topology, no unit filter): 158
  transitions to telomere-like units among 2,253 strongly called species -
  76 one edit, 20 two, 7 three, 55 four or more (18 of these expansions
  containing the old unit). Repeated paths: `AAACCCT->AACCCT` x13 (plants),
  `AACCT->AACCC` x8, `AACCT->ACCTG` x7, `AAACCCT->AACCCTG` x4,
  `AACCT->AAACC` x4. Caveats: topology only, alphabetical tie-breaking (some
  deep inferences and a beetle "reversal" in a polytomy are artefacts);
  Diptera satellite calls to exclude. Next: Mk likelihood with branch
  lengths (B4).
- **G-triplet enrichment** (B3, `repeat_evolution/grun_null.py`). Of 86
  distinct derived telomeric units, 56% have a `GGG` run on the G-rich strand,
  against 26% for the same units shuffled (same length and composition) and
  22% of non-telomeric tandem repeats in the same genomes; `GG` alone is only
  weakly enriched (84% vs 71%). Simulated from the per-unit shuffle
  probabilities: `GGG` p < 1e-5, `GG` p = 1.6e-4. (An earlier "268/268 keep a G-run" was an
  artefact: that count only included units containing `CC`.) To do:
  G-quadruplex potential; telomere-binding protein (TRF/TRB, POT1)
  coevolution.
- **Rates and patterns on the tree.** Number of independent changes,
  reversals, direction, repeated targets (`AACCC`, `ACCTG` in several beetle
  and wasp lineages). Needs a ToL species tree.
- **Mixed telomeres as transition states** (B5,
  `repeat_evolution/mixed_telomeres.py`): 18 species have two telomere-like
  units (<= 12 nt, C-run) at the same chromosome ends, mostly one step apart:
  `AACCT`+`AACCC` (*Bruchidius villosus*, *Oxystoma pomonae*, the ghost moth
  *Korscheltellus lupulinus*), `AACCT`+`AAACC` (*Protonemura*),
  `AACCCACTCT`+`AACCCAGTCT` (*Hedychridium*), `AAACCCT`+`AACCCCT`
  (*Arctium* x2), `AAACCCT`+`AAACCCG` (*Epilobium hirsutum*),
  `AACCCTG`+`AAACCCT` (*Hypericum hirsutum*). Like *Mimulus lewisii* (two TR
  paralogs, mixed telomere; Kumawat 2023). First test (plant TR loci from
  `tr_tert/outputs/tr_template_loci`): ***Hypericum hirsutum* has two TR
  paralogs on different scaffolds, one templating each unit** (SCAFFOLD_5,
  `AAACCCT`; SCAFFOLD_4, `AACCCTG`). *Epilobium hirsutum* (3 copies, one
  cluster) and *Arctium* (4 copies, dispersed) show only `AAACCCT`
  templates. To do: the animal cases (no TRs yet for *Bruchidius*,
  *Oxystoma*, *Korscheltellus*, *Hedychridium*). Caveats: "compound" second units
  (one unit within the other's tandem; 1,280 species) are mostly variant
  repeats inside ordinary telomeres, and many "mixed" second units are
  satellites next to the telomere.
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
- **Tree-wide TERT, first pass** (D1, `tert_tol/`: miniprot with a 505-protein
  panel incl. 159 fungal/protist TERTs, + RBD): found in 2,722, partial 89,
  absent 569 of 3,380 genomes; 319 have neither TERT nor a telomere at any
  end. Diptera 401/402 without TERT (274 also without telomeres); **spiders
  36/36 without TERT** (28 without telomeres) - twice the earlier sample.
  Absences in damselflies (14/14), leafhoppers, mayflies and stoneflies have
  telomeres, so are search misses: profile-HMM + RBD search of the 168
  non-Diptera absences (+20 Diptera controls) running.
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

## Requirements checklist (2026-10-10)

Core = needed for the paper; optional = section 4 or a later paper.

**A. Data and infrastructure**
- [ ] A1 Finish the ToL TR search (round-2/3 stragglers, literature-seeded run) and merge into **one TR table**: species, TR locus, template, copies, model, support
- [ ] A2 TR search for the 7 genomes left out (2 out of memory, 5 too large)
- [ ] A3 Physical consolidation (`tol_telomeres` -> `tr_tert/tol`, `tol_run/`) and path fixes
- [x] A4 **Species tree** (topology): Open Tree of Life induced subtree, 3,363 tips (`tree/species_tree_clean.nwk`); dates still to add - for the 3,422 species: Open Tree of Life synthetic topology + TimeTree dates; BUSCO-based trees for key clades if needed
- [ ] A5 **Final repeat-state table**: canonical unit, support, compound/mixed flag, curated exceptions (*Andrena* dual repeats, satellite calls, *Aphodius*)
- [ ] A6 (optional) Refresh the assembly list; add ERGA/EBP

**B. Section 1, repeat evolution**
- [~] B1 Ancestral-state reconstruction on the tree; transitions counted per lineage, not per species
- [~] B2 Classify transitions: edit distance, "contains ancestral", route 1 vs 2, against a null
- [x] B3 G-run constraint against a null (random units of matched length/composition; non-telomeric tandem repeats from the same genomes)
- [ ] B4 Transition rates (Mk models), repeated targets, reversals
- [x] B5 Mixed/compound telomeres: second unit at chromosome ends, per species
- [ ] B6 Predictors (phylogenetic regression): holocentricity (literature list), genome size, GC, chromosome number, TR copy state

**C. Section 2, TR mechanism**
- [ ] C1 TR table tree-wide (from A1) with template checks
- [ ] C2 Positional survey extended to vertebrates, plants and the new animal groups
- [ ] C3 TR copy number per genome with a uniform search per clade
- [ ] C4 Route-vs-copy test (phylogenetically controlled)
- [ ] C5 TE content around mobile vs stationary TRs (ToL Earl Grey repeat annotations where present)
- [ ] C6 TERT coevolution: codon alignments (CDS from miniprot or ToL gene annotations), PAML branch/clade models
- [ ] C7 Biogenesis classification: type-3 promoter, box H/ACA, intronic

**D. Section 3, telomerase loss**
- [~] D1 **TERT search across all 3,422 species** (miniprot done; HMM fallback running) (miniprot panel + RBD, HMM fallback)
- [ ] D2 Map losses (no TERT, no telomere at ends); characterise chromosome ends without telomerase
- [ ] D3 Spiders: more evidence (sensitivity checks, more assemblies)

**E. Section 4 (optional)**
- [ ] E1 Telomere length atlas (`tidk ends` lengths; `tidk length` on HiFi)
- [ ] E2 Interstitial telomeric sequences as fusion scars
- [ ] E3 TR horizontal-transfer screen
- [ ] E4 T2T QC table for DToL

**F. Figures and writing**
- [ ] F1 Fig 1: tree with repeat states; Fig 2: two routes + G-run; Fig 3: TR validation and relocation; Fig 4: TR copies and mechanism; Fig 5: telomerase loss
- [ ] F2 Methods (tidk 0.3.0 and the TR/TERT pipelines), supplementary tables, data release

Order: start the long compute first (D1, A1), then A4 (needed by B1, B4, B6, C4), then the quick analyses (B2, B3, B5) while it runs.

## Status

Done: 1 (data, first numbers), 2 (validation, relocation survey, copies), 3
(TERT survey). Running: ToL-wide TR search with literature seeds (molluscs,
cnidarians, echinoderms, plants beyond flowering plants, fungi). Missing: a
species tree; per-lineage counts; the tests marked "To do".
