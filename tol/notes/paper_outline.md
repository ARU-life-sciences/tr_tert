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

## Prior work that overlaps, and what stays new (2026-10-11)

- **Fajkus et al. 2026 GBE 18:evag069 (beetles).** TR in *Tribolium* (two
  paralogs, TcTR1/2, active in vitro); a Coleoptera CM (Zenodo
  17782128) finds TRs in 301/420 beetle genomes (Table S2: 299
  sequences), templates matching lineage-specific repeats; **TERT absent
  in derived Scarabaeoidea (Dynastinae, Rutelinae, Cetoniinae, Hopliini)**,
  present in basal Melolonthinae. 184 of our 317 beetles are in their set
  (mostly the same DToL assemblies). Their telomere motifs agree with our
  calls in 159/177 shared species; their TERT calls agree with ours except
  seven rove beetles/histerids where we missed TERT (rerun with their
  beetle TERTs). So: beetle TRs and the scarab loss are theirs - we cite,
  use their CM, and fold beetles into the tree-wide analyses.
- **Lukhtanov & Pazhenkova 2023 Biol J Linn Soc 140:536.** Motif survey of
  220 DToL insects (8 orders): non-canonical motifs in Pentatomomorpha and
  Apocrita; SART/TRAS insertions in telomeres of most insects; long repeats
  in Syrphidae/Tachinidae; hypothesis that a telomerase/retrotransposon
  balance relaxes motif constraint. Our calls agree in 46/51 shared species
  (disagreements are our weak calls). Their hypothesis is testable with our
  data (repeat change rates vs SART/TRAS at ends, TERT loss).
- **Stoianova et al. 2025 Diversity 17:552; 2026 CompCytogen 20:171
  (Hemiptera).** Motifs from 16 assemblies (used tidk) and TREW on SRA
  reads for Nepomorpha: `TTAGAGGTGG` (*Gerris*), `TTAGGGATGG` (Aradidae),
  `TGAC` (*Aphrophora*), 10-bp `TGGATAGG..` decamers in Nepomorpha. No TRs,
  no TERT. Hemiptera TRs remain unknown: open for us.
- **Samo et al. 2026 bioRxiv (*Mimulus*).** TERT evolves fast in lineages
  sharing the TR duplication; TERT binds both TR paralogs. Supports C6 (TERT
  coevolution with TR duplication) as a tree-wide test.
- **Klapproth et al. 2026 (*Caenorhabditis*).** Intronic TR (in *nmy-2*),
  21 species. Another biogenesis switch for C7.

What stays new: the scale and the comparisons across clades - per-lineage
transition counts and rates on a dated tree across all eukaryote groups in
DToL, the two routes and the G-run constraint, TR relocation and copy
number compared between clades, tree-wide telomerase loss (spiders new),
and Hemiptera TRs if found. Motif catalogues per group and single-clade TR
discovery are no longer novel on their own.

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
- [~] A1 Finish the ToL TR search and merge into **one TR table** (`tol_telomeres/shifts/tr_table.tsv`: species, TR locus, template, copies, model, support). 2026-10-12: **1,980 validated TRs** of 2,291 strong calls (tr_tol rounds + literature models 1,818; Coleoptera via the Fajkus et al. 2026 CM 153; round 4 6; Hemiptera 3, Pentatomidae). Template agrees with the call in 2,032/2,057 strongly called species with a template. Two beetles (*Chrysolina hyperici*, *Dascillus cervinus*) still running. Gaps: Hemiptera outside Pentatomidae (3/40), Basidiomycota 0/20, ascidians 0/16, most molluscs/annelids, Odonata, Collembola, crustaceans. Model-free searches failed there: template co-variation within genera (`structcov`; needs a repeat change in the group; *Eurysa*/*Conomelus* leads not TRs, 2026-10-12) and template-bearing RNA-seq reads (`rnatr`; ToL libraries poly-A, TR absent in 6/10 benchmarks). Closing the gaps needs rRNA-depleted RNA-seq. 2026-10-13, Hemiptera: public rRNA-depleted/virome libraries for 6 species (TR reads ~0.04 per million; *Nezara* TR, found by the pentatomid profile, unit `AACCATCCCC` from its template, has 0 template reads in 84 M) and a delphacid conservation + template-switch test (51,146 *Sogatella* loci; Eurysa templates at chance level, 26/11,063) found no TR outside Pentatomidae
- [ ] A2 TR search for the 7 genomes left out (2 out of memory, 5 too large)
- [ ] A3 Physical consolidation (`tol_telomeres` -> `tr_tert/tol`, `tol_run/`) and path fixes
- [x] A4 **Species tree, dated**: Open Tree of Life induced subtree, 3,363 tips; node ages from TimeTree pairwise divergence times (935 of 2,463 internal nodes dated, the rest interpolated BLADJ-style), ultrametric, root 1,516 Myr (`tree/species_tree_dated.nwk`, `tree/date_tree.py`). Checks: Lepidoptera-Hymenoptera 335 Myr, *Bombus*-*Andrena* 99, oak-*Arabidopsis* 107, sawfly-bee 248 (a little young). Caveats: polytomies where the synthesis is unresolved; 62% of node ages interpolated. Was: for the 3,422 species: Open Tree of Life synthetic topology + TimeTree dates; BUSCO-based trees for key clades if needed
- [x] A5 **Final repeat-state table** (`shifts/repeat_evolution/repeat_states.tsv`, 2026-10-12): 3,416 species; 2,678 telomeric, 453 telomerase lost, 285 unknown; 11 calls replaced by the template unit; analysis set 2,664 telomeric (85 units) + 453 lost, 299 excluded (template disagrees, no C run, unknown). Flags: compound, mixed, tert_missed
- [ ] A6 (optional) Refresh the assembly list; add ERGA/EBP

**B. Section 1, repeat evolution**
- [~] B1 Ancestral-state reconstruction on the tree; transitions counted per lineage, not per species
- [~] B2 Classify transitions: edit distance, "contains ancestral", route 1 vs 2, against a null
- [x] B3 G-run constraint against a null (random units of matched length/composition; non-telomeric tandem repeats from the same genomes)
- [~] B4 Transition rates (Mk models), repeated targets, reversals. 2026-10-12: rates done (`b4_rates.R`, `b4_lrt.R`): rate of changing repeat per lineage per 100 Myr, Hymenoptera 0.63, Coleoptera 0.18, Magnoliopsida 0.10, Aves 0.03, Lepidoptera 0.007; Hym/Lep 90-fold (LRT p = 2.5e-36). Repeated targets and reversals still to do (with B1)
- [x] B5 Mixed/compound telomeres: second unit at chromosome ends, per species
- [ ] B6 Predictors (phylogenetic regression): holocentricity (literature list), genome size, GC, chromosome number, TR copy state

**C. Section 2, TR mechanism**
- [x] C1 TR table tree-wide (from A1) with template checks: done with A1 (copy-aware template check per species)
- [x] C2 Positional survey extended to vertebrates, plants and the new animal groups (2026-10-13, `shifts/denovo/outputs/survey_c2_summary.tsv`): 783 new pairs + 250 rerun. TR position kept: birds and mammals 100% to class; fish to family; liverworts, mosses mostly; Lepidoptera to family; flowering plants and Trichoptera within genera only; Hymenoptera 16/25 genus, 0/12 order; Coleoptera 11/25 genus, 1/67 beyond. Moved intervals empty (222/222); another TR copy at the orthologous position in 13 pairs
- [x] C3 TR copy number per genome with a uniform search (`shifts/tr_copies/uniform/copies.tsv`, 2026-10-13): own TR vs own genome, nhmmer E ≤ 1e-10, copy = ≥ 50% coverage. Median 1; 1,344/1,980 single-copy; Lepidoptera mean 1.16 (max 4), Hymenoptera 2.67 (max 69), Coleoptera 3.61 (max 38)
- [~] C4 Route-vs-copy test (phylogenetically controlled). 2026-10-13 (`repeat_evolution/c4_mcmcglmm.R`): MCMCglmm threshold model, derived repeat ~ log2 copies + phylogeny, 1,741 species (182 derived): no copy-number effect (arthropods −0.12, HPD −0.42 to 0.18); phylogenetic h2 0.90-0.95; fewer own-template copies in derived arthropods (−0.56, −1.05 to −0.16). Ancestral-template copies: 55 in 24 derived species (Andrena, Bruchidius, Sitona, Potentilla); counting them removes the templated-copies effect (−0.22, −0.65 to 0.20). Clade level: C2 relocation order = B4 rate order. Still to do: formal clade-level test; route 1 vs 2 per transition
- [ ] C5 TE content around mobile vs stationary TRs (ToL Earl Grey repeat annotations where present)
- [ ] C6 TERT coevolution: codon alignments (CDS from miniprot or ToL gene annotations), PAML branch/clade models
- [ ] C7 Biogenesis classification: type-3 promoter, box H/ACA, intronic

**D. Section 3, telomerase loss**
- [x] D1 **TERT search across all 3,422 species**: found 2,749, found by HMM fallback 80 (Odonata 14, Hemiptera 12), partial 91, absent 492, not run 4 (very large genomes); absent = Diptera 402/402 (the *Tipula helvola* hit is a microsporidian contaminant, *Nosema bombycis*), spiders 38/38, ~53 scattered others (Staphylinidae, leafhoppers, nemerteans, *Littorina*), most with canonical telomeres, so probably divergent TERTs missed (miniprot panel + RBD, HMM fallback)
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

2026-10-13: data collection essentially done (A1, A4, C1, D1). Done: A5, B3, B5, C2, C3; partial: B1, B2, B4
(rates done), C4 (no species-level copy effect; lineage-level match of TR relocation and repeat-change rate).
Hemiptera TR search outside Pentatomidae closed as a negative result.
Next: a formal clade-level test of rate vs relocation (C4), D2, B6, B1/B2 on the dated tree. Prose draft:
`paper_draft.md` (`[TODO]` markers point to these items).
