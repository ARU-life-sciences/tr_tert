# Draft: How telomeres evolve across the Tree of Life (2026-10-12)

First prose draft of the phase 4 paper, built on `paper_outline.md` (plan,
checklist) and `literature_synthesis.md` (prior work). Numbers are current
as of 2026-10-12; `[TODO ...]` marks results still to come, with the
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
called repeat in 98.8% of species, and TRs duplicate and move frequently in
some lineages (Hymenoptera: 3.3 TR-like copies per genome, positions not
conserved between families) but rarely in others (Lepidoptera: 1.2 copies,
positions conserved). [TODO C4: whether large changes coincide with
duplicated or relocated TRs.] Telomerase was lost in all Diptera and all
spiders (38/38, new), where chromosome ends have no telomeric repeat.
[TODO: one-sentence conclusion once C4 and B4 are in.]

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
  Hymenoptera (40) and spiders (30). [TODO A5: final repeat-state table
  with curated exceptions; Fig 1 tree.]
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
  [TODO B4: Mk rates with branch lengths on the dated tree; reversals;
  per-lineage rates. TODO B2: null for the edit-distance distribution.
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
  polyadenylated, is missing from 6 of 10 benchmark species). [TODO:
  planthopper lead, `structcov/leads/`; rRNA-depleted libraries would
  close the gaps.]
- **TRs duplicate.** TR-like copies (E ≤ 1e-5): Hymenoptera mean 3.3 per
  genome, Lepidoptera 1.2 (p = 3e-22); pimpline wasps 3-7 copies with one
  carrying a template; halictid bees up to 52. [TODO C3: uniform copy
  search per clade.]
- **TRs move.** In 250 congeneric and confamilial pairs, the TR is at the
  orthologous position in Lepidoptera (genus 25/25, family 22/25) but often
  not in Hymenoptera (16/25, 4/24); where it moved, the orthologous interval
  is conserved and empty, as in *Mimulus*. [TODO C2: vertebrates, plants,
  beetles.]
- **How TRs move.** No retrocopy signature (poly-A tails and target-site
  duplications at background rates). Copies are either dispersed TR-sized
  insertions with unique flanks (pimplines, *Lasioglossum*) or copies
  inside a larger repeated unit (*Andrena*, *Seladonia*). [TODO C5: TE
  annotation of flanks, copy ages.]
- **Do mobile TRs make big changes?** [TODO C4, the key test: large
  repeat changes (route 2) on branches with multi-copy or relocated TRs,
  small ones (route 1) with single-copy TRs; phylogenetically controlled.]
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
2. TR mobility differs by orders of magnitude between lineages; if C4
   holds, TR duplication and relocation explain why some lineages
   (Hymenoptera, Hemiptera?) have diverse repeats and others (Lepidoptera)
   do not. This generalises *Mimulus* and *Andrena* to the tree.
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
- **Repeat transitions.** Fitch parsimony; edit distance between parent
  and child units; containment of the ancestral unit. [Mk models, B4.]
- **G-run null.** Each derived unit shuffled 1,000 times (same length and
  composition); non-telomeric tandem repeats from the same genomes.
- **TR search.** nhmmer family profiles from published TRs, Infernal CMs
  (Rfam RF00024; sawfly, beetle CMs), model walking with validated seeds;
  copy-aware template check (exact match ≥ unit + 2 nt to the tandem, in
  any copy); validation rule as above.
- **TR copies and position.** TR-like loci at E ≤ 1e-5, one per 2 kb;
  synteny of TR flanks between pairs of species; poly-A, TSD and TATA
  tests against random positions.
- **TERT.** miniprot (v0.18), panel grown over rounds from 40 UniProt TERTs
  to 505 proteins; Telomerase_RBD (PF12009) hmmsearch; profile HMM from 236
  TERTs on six-frame ORFs, positive at RBD i-E ≤ 1e-4.
- **Searches that did not find TRs** (supplement): template co-variation
  and RNAz/R-scape within genera (`structcov`); template-bearing RNA-seq
  reads (`rnatr`) with ORF and protein triage calibrated on 60 known TRs.

## Figures

1. Dated tree of 3,415 species: repeat state, TR found, TERT state.
2. Two routes: edit-distance distribution of transitions; recurrent edits;
   G-triplet enrichment against the nulls.
3. TR validation (template vs call) and TR position conservation
   (Lepidoptera vs Hymenoptera).
4. TR copies by clade; mixed telomeres and *Hypericum* paralogs; the C4
   test.
5. Telomerase loss on the tree; chromosome ends without telomerase.
