# Draft: How telomeres evolve across the Tree of Life (2026-10-13)

First prose draft of the phase 4 paper, built on `paper_outline.md` (plan,
checklist) and `literature_synthesis.md` (prior work). Numbers are current
as of 2026-10-13; `[TODO ...]` marks results still to come, with the
checklist item from the outline. Working title:

*How telomeres evolve: repeat changes, telomerase RNA mobility and
telomerase loss across 3,400 genomes*

## Abstract (provisional)

Telomeric repeats are among the most conserved sequences in eukaryote
genomes, yet many lineages have changed theirs, and some have lost
telomerase altogether. Because each change must be written by the
template of the telomerase RNA (TR), repeat evolution is a record of TR
evolution, but TRs are short, fast-evolving and hard to find, and have
been studied one genus or phylum at a time. We called telomeric repeats at
the chromosome ends of 3,415 chromosome-level genomes from the Darwin Tree
of Life, found and validated TRs in 1,980 of them, and searched all for the
telomerase reverse transcriptase (TERT). Repeats change by two routes: [N]
single template edits, and [N] larger changes that usually keep the
ancestral repeat and add to it. Derived repeats are strongly enriched for
G-triplets (56% vs 26% of shuffled units). The TR template encodes the
called repeat in 98.8% of species. The rate of repeat change differs
90-fold between orders on a dated tree (Hymenoptera 0.63 changes per
lineage per 100 Myr, Lepidoptera 0.007), and so do TR duplication and
movement: the TR has kept its genomic position across all birds and
mammals and across Lepidoptera families, but moved in most beetle and
Hymenoptera pairs beyond the genus, with the old position left empty. Within
clades, however, a species' TR copy number does not predict whether its
repeat has changed: derived repeats cluster in a few lineages
(phylogenetic heritability 0.9). Telomerase was lost in all Diptera and
all spiders (38/38, new), where chromosome ends have no telomeric repeat.
Lineages whose TR moves are the lineages whose telomeric repeat changes.

## Introduction

1. **Telomeres and telomerase.** Linear chromosomes end in tandem repeats
   synthesised by telomerase, a reverse transcriptase (TERT) that copies a
   short template in its RNA subunit (TR). The repeat is therefore set by
   the template, and a change of repeat requires a change of template.
2. **The repeat is conserved, but not everywhere.** `TTAGGG` in vertebrates
   and many other animals, `TTTAGGG` in most plants, `TTAGG` in arthropods;
   but many lineages have shifted (beetles, Hymenoptera, Hemiptera, some
   plants), and some have abandoned telomerase (Diptera use
   retrotransposons or satellites). Previous surveys catalogued motifs in
   tens to a few hundred genomes per group (Lukhtanov & Pazhenkova 2023,
   220 insects; Stoianova et al. 2025, 2026, Hemiptera; Fajkus et al. 2026,
   beetles), without TRs or with TRs in one clade.
3. **TRs are hard to find**, so the mechanism of repeat change has been
   studied in single genera. TR discovery has gone lineage by lineage, by
   iterative covariance-model searches outward from a known TR
   (vertebrates, echinoderms, eight more animal phyla, plants, Hymenoptera,
   *Caenorhabditis*, beetles). In *Mimulus* (Kumawat et al. 2023) and
   *Andrena* (Klapproth et al. 2025), TR paralogs with different templates
   make sequence-heterogeneous telomeres, and the TR moves between genomic
   positions; Kumawat et al. proposed "transposition, duplication and
   divergence" as a general model.
4. **What is unknown**: how often repeats change, by what kinds of
   mutation, whether some lineages change more because their TRs duplicate
   and move, and where telomerase has been lost. Answering this needs
   repeats, TRs and TERT across the same set of genomes, on a dated tree.
5. **This study.** The Darwin Tree of Life has produced chromosome-level
   assemblies for thousands of eukaryotes, from the same sequencing and
   curation pipeline. We call telomeres at chromosome ends with tidk 0.3.0,
   find TRs by model walking with telomere-validated seeds, and search
   TERT, for 3,415 species across 19 phyla.

## Results

### 1. Telomeric repeats across 3,415 genomes

- **Calls.** tidk 0.3.0 (error-tolerant runs, primitive canonical units,
  chromosome-end checks) gives a strong telomere call (the unit at a
  substantial fraction of chromosome ends) in 2,291 species, weak in 692,
  none in 432. Species without a call are mostly Diptera (274), fungi (52),
  Hymenoptera (40) and spiders (30).
- **Repeat states.** Combining calls, TR templates and TERT
  (`repeat_states.tsv`): 2,678 species telomeric, 453 telomerase lost
  (TERT absent, with no strong call, or in Diptera and spiders), 285
  unknown. Where a strong call disagrees with the TR template and the
  template unit is at ≥ 4 chromosome ends, the template unit is taken (11
  species, e.g. *Sabra* `ACCTC` → `AACCT`); 9 other disagreements and
  units with no C run and no template are excluded from analyses. Analysis
  set: 2,664 telomeric species with 85 distinct units, 453 losses. [Fig 1
  tree.]
- **Calls are confirmed by the TR template** (section 2): where a TR is
  found, its template encodes the called unit in 2,032 of 2,057 strongly
  called species (98.8%). Most of the 21 disagreements are tidk miscalls
  (a satellite or a unit one step off), not template disagreements.
  Independent checks: our calls agree with Lukhtanov & Pazhenkova (2023)
  in 46/51 shared species and with Fajkus et al. (2026) in 159/177 beetles,
  the disagreements being our weak calls.
- **Two routes of change.** By Fitch parsimony on the Open Tree of Life
  topology, 158 transitions to telomere-like units among 2,253 strongly
  called species: 76 are one edit from the ancestral unit, 20 two, 7 three,
  and 55 four or more, with few in between. Most large changes keep the
  ancestral unit and add to it (`AACCT` → `AACCCAGACCT`, `AACCATCCCT`;
  18 of 55 contain the old unit exactly). The same edits recur:
  `AAACCCT` → `AACCCT` 13 times in plants, `AACCT` → `AACCC` 8 times,
  `AACCT` → `ACCTG` 7 times, `AACCT` → `AAACC` 4 times.
- **Rates differ by orders of magnitude.** An equal-rates Mk model of the
  repeat state on the dated tree, per clade, gives the rate at which a
  lineage changes repeat: Hymenoptera 0.63 per 100 Myr (95% CI
  0.50-0.80), Coleoptera 0.18 (0.13-0.24), Hemiptera 0.10, flowering
  plants 0.10 (0.06-0.16), birds 0.03, Lepidoptera 0.007 (0.001-0.022).
  Hymenoptera vs Lepidoptera 90-fold (likelihood ratio test of a shared
  rate, p = 2.5e-36); Coleoptera vs Lepidoptera 26-fold (p = 6e-13);
  Hymenoptera vs Coleoptera 3.5-fold (p = 9e-11). Parsimony steps per
  1,000 Myr of tree give the same order (Hymenoptera 5.3, Coleoptera 1.7,
  Lepidoptera 0.07). [TODO B1/B2: reversals and repeated targets on the
  dated tree; null for the edit-distance distribution.
  Caveat to fix: alphabetical tie-breaking in polytomies, Diptera
  satellite calls.]
- **G-triplet constraint.** Of 86 distinct derived units, 56% have a `GGG`
  run on the G-rich strand, against 26% of the same units shuffled (same
  length and composition; p < 1e-5) and 22% of non-telomeric tandem
  repeats from the same genomes. `GG` alone is only weakly enriched (84% vs
  71%, p = 1.6e-4). New repeats keep the feature needed for G-quadruplex
  formation or telomere-protein binding. [TODO: G4 potential; TRF/POT1.]
- **Mixed telomeres as intermediate states.** 18 species have two
  telomere-like units at the same chromosome ends, mostly one edit apart
  (`AACCT` + `AACCC` in two beetles and a ghost moth; `AAACCCT` +
  `AACCCCT` in two *Arctium*). In *Hypericum hirsutum* the two units are
  templated by two TR paralogs on different scaffolds (`AAACCCT` and
  `AACCCTG`), as in *Mimulus lewisii*: a duplicated TR whose copies
  diverged in template can make a mixed telomere, a possible step towards a
  new repeat.
- [TODO B6: predictors - holocentricity, genome size, GC, karyotype, TR
  copy state - by phylogenetic regression.]

### 2. Telomerase RNA: validation, duplication and movement

- **A tree-wide TR table.** Starting from published TRs and Rfam models
  (Hymenoptera, Lepidoptera, vertebrates, plants, beetles; Fajkus et al.
  2021, 2023, 2026; Xie et al. 2008; Logeswaran et al. 2021), we searched
  every genome with lineage-specific profiles and covariance models,
  walking outward: each new TR whose template matched its species'
  telomere became a seed for the next round. A TR is validated when the
  species has a strong telomere call, the template matches it (at least
  the unit length + 2 nt of its tandem), and the hit has E ≤ 1e-10:
  **1,980 validated TRs** (Lepidoptera 779, flowering plants 476,
  Coleoptera 154, Hymenoptera 155, birds 110, mosses 109, fish 35,
  mammals 26, liverworts 30, plus molluscs, echinoderms, annelids, green
  algae, ...). [Table 1: by clade.]
- **Gaps.** No validated TR in Hemiptera outside Pentatomidae (3/40),
  Basidiomycota (0/20), ascidians (0/16), most molluscs and annelids,
  Odonata, Collembola, crustaceans. These groups have telomeres and TERT
  but no seed close enough for model search. Two searches without a model
  failed: template co-variation within genera (works only where a genus
  includes a repeat change, as in *Bombus*), and template-bearing reads in
  RNA-seq (the ToL libraries are poly-A selected, and the TR, which is not
  polyadenylated, is missing from 6 of 10 benchmark species). Public
  rRNA-depleted and virome libraries do not help: TR reads are about 0.04
  per million (3 in 72 M in *Scotinophara*), and the *Nezara viridula* TR
  (found by the pentatomid profile; its template predicts `AACCATCCCC`)
  has no template read in 84 M. A conservation test in Delphacidae (51,146
  transcribed template loci of *Sogatella*, homologs in *Laodelphax*,
  *Nilaparvata*, *Conomelus* and *Eurysa*, which has a different repeat)
  found no locus whose template follows the repeat change: Eurysa
  homologs carry their own template no more often than by chance (26 of
  11,063). Pentatomid profiles, relaxed covariance models and Eurysa
  co-variation leads found nothing outside Pentatomidae.
- **TRs duplicate.** One uniform search per species (its own TR as the
  query, E ≤ 1e-10, copies = hits covering ≥ half the TR): most TRs are
  single-copy (1,344 of 1,980; median 1), but Hymenoptera average 2.7
  copies (up to 69, *Tiphia femorata*) and Coleoptera 3.6 (up to 38,
  *Amara*), against 1.2 in Lepidoptera (max 4). Most extra copies carry a
  template. Pimpline wasps have 3-7 copies with one carrying a template;
  halictid bees up to 47.
- **TRs move, in some lineages.** For 1,033 pairs of species with
  validated TRs (up to 25 per clade and taxonomic level), the interval
  orthologous to species A's TR was projected into species B by
  nucleotide and BUSCO-gene anchors. The TR stays at the orthologous
  position in every bird and mammal pair, up to different orders (25/25
  and 25/25 at class level), in fish to family (10/10), in liverworts and
  mostly in mosses (23/25 at order level), and in Lepidoptera to family
  (22/25). It is stable within genera but mobile beyond in flowering
  plants (genus 23/24, subfamily 7/16, family 3/11) and caddisflies (genus
  5/5, order 4/13). It moves often in Hymenoptera (genus 16/25, family
  4/24, order 0/12) and most of all in beetles (genus 11/25, subfamily
  1/18, family and beyond 0/49). Where the TR moved, the orthologous
  interval holds no TR-like sequence (no hit to B's TR at E ≤ 1e-10 in any
  of 222 moved pairs): relocation, not decay in place, as in *Mimulus*.
  Another copy of B's TR sits at the orthologous position in only 13
  pairs.
- **How TRs move.** No retrocopy signature (poly-A tails and target-site
  duplications at background rates). Copies are either dispersed TR-sized
  insertions with unique flanks (pimplines, *Lasioglossum*) or copies
  inside a larger repeated unit (*Andrena*, *Seladonia*). [TODO C5: TE
  annotation of flanks, copy ages.]
- **Do species with more TRs change repeat?** Between clades, yes: the
  orders with the fastest repeat change (Hymenoptera, Coleoptera) have the
  most TR copies. Within clades, no: a phylogenetic probit model
  (MCMCglmm, dated tree as random effect; 1,741 arthropods and plants, 182
  with a derived repeat) finds no effect of TR copy number (log2 copies,
  arthropods −0.12, 95% HPD −0.42 to 0.18; Hymenoptera + Coleoptera −0.11,
  −0.43 to 0.19). Derived repeats cluster by lineage (phylogenetic
  heritability 0.90-0.95), so the association is at the level of clades,
  not species. An apparent deficit of templated copies in derived
  arthropods (−0.56, −1.05 to −0.16) disappears when copies carrying the
  *ancestral* template are counted (−0.22, −0.65 to 0.20): 55 copies in 24
  derived species kept only the old template, most in *Andrena* (up to 11
  in *A. bucephala*), *Bruchidius*, *Sitona* and *Potentilla*: old TR
  copies left behind after the repeat changed.
- **Mobility and repeat change go together across lineages.** Ranking
  clades by TR relocation (C2) gives the same order as their rate of
  repeat change (B4): beetles and Hymenoptera highest, flowering plants
  and caddisflies intermediate, Lepidoptera, birds and mammals lowest.
  [TODO: formal clade-level test (rate vs fraction relocated); route 1 vs
  2 per transition.]
- [TODO C6: TERT rates on branches with TR duplication or relocation.
  TODO C7: biogenesis class (type-3 promoter, box H/ACA, intronic) mapped
  on the tree.]

### 3. Telomerase loss

- **TERT across the tree.** miniprot with a 505-protein panel (including
  159 fungal and protist TERTs), confirmed by the TERT RNA-binding domain,
  with a profile-HMM fallback: TERT found in 2,832 genomes (80 by the HMM
  only), partial in 91, absent in 488.
- **Losses with no telomere.** 314 species have neither TERT nor a
  telomeric repeat at any chromosome end. All 402 Diptera lack TERT (274
  also without telomeres; one apparent TERT, in *Tipula helvola*, is a
  microsporidian contaminant). **All 38 spiders lack TERT, 30 without
  telomeres at any end**, extending the suggestion of Fajkus & Fajkus
  (2025) to an order-wide loss. Beetles: the scarab losses of Fajkus et
  al. (2026). Scattered absences in Hemiptera (8/72), Coleoptera (5),
  Hymenoptera (5), nemerteans and *Littorina* are in species with
  canonical telomeres, so are probably divergent TERTs missed by the
  search.
- [TODO D2: map losses on the dated tree; what occupies chromosome ends
  without telomerase (retrotransposon arrays, satellites). TODO D3:
  spider sensitivity tests.]
- **Shift lineages keep telomerase.** In 30 lineages with a changed
  repeat, 196/198 species have TERT: repeat changes are telomerase
  products, not a switch to a telomerase-independent mechanism.

### 4. Further results (supplement or separate)

Telomere length atlas from `tidk ends` and `tidk length`; interstitial
telomeric sequences as fusion scars; a weevil TR 96% identical to a
riodinid butterfly TR (horizontal transfer?); T2T status of every DToL
assembly. [TODO E1-E4, optional.]

## Discussion (points to make)

1. Telomeric repeats change by two routes, with different mechanisms: a
   template point mutation, and a template expansion that keeps the old
   unit. The G-triplet constraint shapes both.
2. Repeat-change rates and TR mobility both differ by orders of magnitude
   between lineages, in the same direction (Hymenoptera, Coleoptera high;
   Lepidoptera low). The link is a property of lineages, not of
   individual species' copy number: a clade with mobile TRs changes
   repeat more often, but the species that changed need not carry more
   copies today. This generalises *Mimulus* and *Andrena* to the tree.
3. Telomerase loss is rare but recurrent (Diptera, spiders, scarabs), and
   always with loss of the telomeric repeat at the ends.
4. Limits: TR discovery is still seed-bound, so the TR table is biased to
   clades with known TRs; Hemiptera, molluscs, annelids and fungi need
   rRNA-depleted RNA-seq. Parsimony on a partly unresolved, partly
   interpolated tree.
5. Resources: the repeat table, the TR table (1,980 validated TRs, the
   largest set to date), TERT calls, and tidk 0.3.0.

## Methods

- **Genomes.** 3,415 chromosome-level DToL assemblies (curated primary or
  haplotype 1), [date of freeze].
- **Telomere calls.** tidk 0.3.0 `explore` (error-tolerant runs seeded on
  exact adjacent copies; primitive canonical units collapsing rotations,
  reverse complements and multimers; strand balance), on chromosome ends,
  and `ends` for per-end presence and T2T status. Simulations (0-10%
  substitution and indel error): true unit top in 21/21 conditions.
  Strong/weak support thresholds: [define].
- **Species tree.** Open Tree of Life induced subtree (3,363 tips); node
  ages from TimeTree pairwise divergences (935 of 2,463 internal nodes
  dated, the rest interpolated), ultrametric, root 1,516 Myr.
- **Repeat states.** Rules in order: telomerase lost (TERT absent in
  Diptera/Araneae, or absent with no strong call); curated exceptions;
  strong call with agreeing template; strong call replaced by the template
  unit when it is at ≥ 4 ends; weak call with agreeing template; else
  unknown (`repeat_evolution/repeat_states.py`).
- **Repeat transitions.** Fitch parsimony; edit distance between parent
  and child units; containment of the ancestral unit. Rates: phytools
  `fitMk` equal-rates model per clade on the dated tree (polytomies
  resolved at random, 0.01 Myr branches); rate of change r = (k − 1) q;
  likelihood ratio test of one shared r between two clades; profile 95%
  CIs (`b4_rates.R`, `b4_lrt.R`).
- **TR copies vs repeat change.** MCMCglmm threshold model, derived repeat
  (differs from the phylum's ancestral unit) ~ log2 TR copies, phylogenetic
  random effect from the dated tree (inverseA), residual variance fixed at
  1, parameter-expanded prior, 200,000 iterations (`c4_mcmcglmm.R`).
- **G-run null.** Each derived unit shuffled 1,000 times (same length and
  composition); non-telomeric tandem repeats from the same genomes.
- **TR search.** nhmmer family profiles from published TRs, Infernal CMs
  (Rfam RF00024; sawfly, beetle CMs), model walking with validated seeds;
  copy-aware template check (exact match ≥ unit + 2 nt to the tandem, in
  any copy); validation rule as above.
- **TR copies and position.** Per species, its validated TR locus searched
  against its genome (nhmmer, E ≤ 1e-10), hits within 2 kb merged; a copy
  covers ≥ 50% of the query; template-bearing if the hit ± 20 nt matches
  ≥ unit + 2 nt of the tandem (`tr_copies/src/uniform_copies.py`); synteny of TR flanks between pairs of species; poly-A, TSD and TATA
  tests against random positions.
- **TERT.** miniprot (v0.18), panel grown over rounds from 40 UniProt TERTs
  to 505 proteins; Telomerase_RBD (PF12009) hmmsearch; profile HMM from 236
  TERTs on six-frame ORFs, positive at RBD i-E ≤ 1e-4.
- **Searches that did not find TRs** (supplement): template co-variation
  and RNAz/R-scape within genera (`structcov`); template-bearing RNA-seq
  reads (`rnatr`) with ORF and protein triage calibrated on 60 known TRs;
  public rRNA-depleted libraries for 6 hemipterans (genome units from WGS
  reads where assemblies lack telomeres); the delphacid conservation and
  template-switch test (`rnatr/src/planthopper_conservation.py`).

## Figures

1. Dated tree of 3,415 species: repeat state, TR found, TERT state.
2. Two routes: edit-distance distribution of transitions; recurrent edits;
   G-triplet enrichment against the nulls.
3. TR validation (template vs call) and TR position conservation
   (Lepidoptera vs Hymenoptera).
4. TR copies by clade; mixed telomeres and *Hypericum* paralogs; the C4
   test.
5. Telomerase loss on the tree; chromosome ends without telomerase.
