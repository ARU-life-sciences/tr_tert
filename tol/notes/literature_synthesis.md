# Telomerase RNA literature: synthesis and open questions (2026-10-10)

From the 15 papers in `paperpile-files.zip` (repository root), read for what
bears on the ToL dataset.

## What the literature establishes

**TR discovery is lineage by lineage.** Each group's TR was found separately,
mostly by iterative covariance-model searches outward from a known TR:
vertebrates (Xie et al. 2008 x2), echinoderms (Li et al. 2013; Podlevsky et
al. 2016), 82 TRs from 8 more animal phyla (Logeswaran et al. 2021), plants
and other Diaphoretickes (Fajkus et al. 2021), Hymenoptera (Fajkus et al.
2023), Caenorhabditis (Klapproth et al. 2026). None in arthropods outside
Hymenoptera and Lepidoptera, nor in flatworms (Logeswaran et al. 2021).

**Structure is conserved where sequence is not.** A template-pseudoknot
domain and a template boundary element are near universal; a second
TERT-binding domain (CR4/5, or eCR4/5 in echinoderms and protostomes) arose
early (Podlevsky et al. 2016, trypanosome; Podlevsky & Chen 2016).

**Biogenesis switched repeatedly.** Animal TRs are monophyletic and
snoRNA-like (box H/ACA; Logeswaran et al. 2021); plant and ciliate TRs are Pol
III snRNA-like with a type-3 promoter (Fajkus et al. 2021); Hymenoptera
switched to Pol III; *Caenorhabditis* TR sits in an intron of *nmy-2*, from a
single origin (Klapproth et al. 2026; Takeda et al. 2025); fungal TRs are
mRNA-like (Smekalova et al. 2012).

**TRs duplicate and move.**
- *Mimulus* (Kumawat et al. 2023): TR syntenic within sections but not
  between them, with the orthologous interval conserved and empty; paralogs
  with different templates give a sequence-heterogeneous telomere
  (AAACCG + AAACCCG); DNA transposons within 2 kb of the TR in 6/9 species;
  proposed "transposition, duplication and divergence" model of telomere
  evolution.
- *Andrena* (Klapproth et al. 2025): 2+ TR paralogs in most species, in
  syntenic clusters, two templates co-expressed (5' RACE), making two repeats.
- Plants generally: TR paralogs (Závodník et al. 2023, cited).

**TERT-TR coevolution is unresolved.** In *Mimulus*, TERT evolves rapidly
(PAML clade model) in lineages sharing the TR duplication and binds both
paralogs, and *M. lewisii* telomerase is dominant in hybrids (Samo et al.
2026); in *Andrena*, TERT is single copy with no corresponding change
(Klapproth et al. 2025).

**TR size tracks genome size** in fish (teleost TRs 312-348 nt; R² = 0.50
against genome size; Xie et al. 2008).

**TRs can cross genomes.** Marek's disease virus carries a functional TR 88%
identical to chicken TR (Fragnet et al. 2005).

**Telomerase can be dispensed with.** ALT (homology-directed telomere
synthesis, lagging-strand bias) maintains telomeres without telomerase in
~15% of cancers (Min et al. 2017).

## What our dataset adds

Every study above works on one genus or one phylum, with tens of genomes. We
have, for ~3,400 DToL species with chromosome-level assemblies: telomere calls
checked at chromosome ends; TERT in ~330 species (shift lineages, controls,
telomerase-loss tests); validated TRs in ~1,600 species (Lepidoptera,
Hymenoptera, vertebrates, flowering plants; more groups from the literature
seeds now running); TR copy number and position for every one; a 250-pair
positional survey.

## Questions we can now answer

1. **How often do TRs move, tree-wide, and why only in some lineages?**
   *Mimulus* and *Andrena* each show relocation in one genus; we see it
   across Hymenoptera but not Lepidoptera (genus 16/25 vs 25/25 conserved).
   Extend the survey to vertebrates, plants, and the new animal groups, and
   ask what predicts mobility: karyotype stability, genome size, TE load,
   biogenesis class (Pol III snRNA-like vs H/ACA snoRNA-like).
2. **What moves them?** Kumawat proposed DNA transposons (6/9 *Mimulus* TRs
   have one within 2 kb). Test at scale: TE content around mobile vs
   stationary TRs. Our first look found no retrocopy signature and two copy
   patterns (dispersed TR-sized insertions; copies inside larger repeats).
3. **Do TR duplications drive telomeric-repeat shifts?** The hypothesis of
   both Kumawat and Klapproth, so far argued from one genus each. With ~30
   independent shift lineages, matched controls and TR copy numbers, test
   whether shifts are concentrated where TRs are multi-copy or have moved,
   and whether species with two functional templates have mixed telomeres
   (tidk compound calls: *Andrena*, *Carlina*, *Bombus*).
4. **Does TERT coevolve with TR changes?** *Mimulus* says yes, *Andrena* no.
   We have TERT proteins for ~300 shift and control species: dN/dS clade
   models on branches with TR duplication, relocation or a repeat shift.
5. **Does TR size track genome size beyond fish?** Test across vertebrates,
   Lepidoptera, Hymenoptera and plants (needs TR boundaries: CM-aligned
   lengths).
6. **Where did TR biogenesis switch?** Classify each validated TR (type-3
   promoter, box H/ACA, intronic) and map the switches on the tree; look for
   *Caenorhabditis*-like intronic TRs elsewhere.
7. **How often is telomerase lost, and what replaces it?** TERT and TR absence
   plus no telomere at chromosome ends (Diptera, spiders, *Viscum*, ...):
   map the losses, and characterise what is at chromosome ends instead
   (retrotransposon arrays, satellites, ALT-like signatures).
8. **Do TRs move between genomes?** The weevil "TR" 96% identical to a
   riodinid butterfly's, and the viral precedent: screen all genomes for TRs
   whose best match is a distant clade.

Strongest for a paper, in our view: 1 + 3 (relocation and duplication as the
engine of telomere evolution, tested tree-wide rather than in one genus),
with 4 and 7 as major supporting analyses.
