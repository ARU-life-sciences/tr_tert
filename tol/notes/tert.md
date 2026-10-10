# TERT in the telomeric-repeat shift lineages (2026-10-09)

Do the animal lineages whose telomeric repeat changed still have
telomerase? For the 30 shift lineages from the ToL-wide tidk run (a
non-ancestral repeat at >= 25% of chromosome ends, in >= 2 species of an
order), up to 5 controls each from the same family or order with the
ancestral repeat, and three test groups: Diptera (telomerase lost),
spiders (no telomere at any chromosome end) and mayflies (no telomere
call). 329 species; `src/select_species.py` -> `inputs/species.tsv`.

## Method

1. **miniprot** (v0.18) spliced alignment of a TERT protein panel, with a
   Telomerase_RBD (Pfam PF12009) hmmsearch of each translated alignment.
   The RBD is what tells TERT from retrotransposon reverse transcriptases,
   which share its RT domain. Three rounds, the panel growing each time:
   40 UniProt TERTs (`inputs/tert_panel.tsv`), then + the 166 intact DToL
   TERTs from round 1, then + 192 from round 2 and four more UniProt
   proteins (a planthopper, two *Caenorhabditis*, a spider fragment).
   Not-found species fell from 115 to 34 over the rounds, mostly clades
   with no panel relative at first (Carabidae 18/19 not found in round 1).
2. **Profile HMM** for the species still not found: an HMM built from all
   236 panel + DToL TERTs (MAFFT, hmmbuild; 788 match states), and the RBD
   HMM, against all six-frame ORFs >= 40 codons. The summed TERT-profile
   score does not discriminate (Diptera reach 2,235 from retrotransposon RT
   families); an RBD hit at i-Evalue <= 1e-4 does: 0/6 Diptera, 6/8
   miniprot-positive species.

Scripts: `src/run_tert.bash` (miniprot), `src/expand_panel.py`,
`src/summarise_tert.py`, `src/build_tert_hmm.bash`, `src/run_tert_hmm.bash`,
`src/summarise_tert_hmm.py`, `src/combine_tert.py` ->
`outputs/tert_calls.tsv`, `outputs/tert_lineages.tsv`.

## Result

| group | TERT found |
|---|---|
| shift species (30 lineages) | **196/198** |
| controls (ancestral repeat) | 105/108 |
| Diptera | **0/6** |
| spiders | **0/6** |
| mayflies | 5/6 |

- **Every shift lineage keeps TERT.** The repeat changes are telomerase
  products, not a switch to a telomerase-independent mechanism, so they
  are a test of TR template change. The two shift species not found are
  *Chrysis angustula* (2/3 of its lineage found) and *Caenorhabditis
  monodelphis* (whose AAGCCT/TTAGGC is the ancestral nematode repeat, not
  a shift: the selection's ancestral rule of AACCCT for non-arthropod
  animals is wrong for nematodes).
- The 3 controls not found are leafhoppers (Cicadellidae), which have
  canonical telomeres, so they are search misses: 5/8 leafhoppers were
  only found by the HMM. Sensitivity is lowest for clades far from every
  panel protein.
- **Spiders: no TERT by any search**, like Diptera, and no telomere at any
  chromosome end. Consistent with telomerase loss in Araneae (suggested
  in Fajkus & Fajkus 2025). Caveat: the HMM finds 6/8 known positives, so
  a very divergent spider TERT could be missed; the only spider TERT in
  UniProt is short fragments (*Caerostris*, <= 363 aa), and its use as a
  query found nothing.
- "found_fs" (100 species: the best alignment has 1-3 frameshifts or stops
  over >= 90% of a reference) is mostly alignment or assembly error, not
  pseudogenisation: 35 controls with working telomeres are in this class.

## Next

- TR: Hymenoptera TRs are known (Fajkus et al. 2023), so the six
  Hymenoptera shift lineages are where a TR template can be checked against
  the repeat directly. Elsewhere TR needs the de novo pipeline
  (`../../tr_tert/clade_discovery`).
- TERT sequence: shift-associated changes (e.g. in the TRBD or motif T)
  across the 196 shift and 105 control proteins (`outputs/tert_round*_proteins.faa`).
- Spiders: more spider assemblies, and a spider-specific search (e.g.
  tblastn of the RT motifs only).
