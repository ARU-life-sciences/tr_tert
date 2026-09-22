# Clade discovery: extending TR/TERT identification beyond plants

## Motivation

The parent project (`../`) identifies telomerase RNA (TR), TERT, and
telomeric repeats across ~700 DToL plant genomes. Two papers shared by
the user reframe why this might generalize far beyond plants:

- **Fajkus et al. 2023, NAR** (10.1093/nar/gkac1202) - found that
  Hymenoptera (and candidate Lepidoptera) TRs are *not* the expected
  Animalia-type H/ACA-box, RNAPII-transcribed snoRNA - they're short,
  RNAPIII-transcribed, type-3-promoter snRNAs, the *same* biogenesis
  class as plants and ciliates. An independent evolutionary switch,
  not inheritance - but it means the plant-side toolkit/mindset (template-
  region-as-anchor, telomere-motif-variability-as-signal) has a real
  chance of transferring to Arthropoda specifically, unlike "typical"
  Animalia.
- **Fajkus & Fajkus 2025, Biochem J review** (10.1042/BCJ20240501) -
  explicitly names the priority unknowns: "Understanding unknown TR
  characteristics in Ecdysozoa (e.g. nematodes and arthropods) would
  substantially enhance our knowledge..." Also names specific
  telomerase-loss candidates worth checking: Cynipoidea (already
  confirmed lost in the 2023 paper), spiders (Araneae, citation [43]),
  *Meloidogyne* (Nematoda), and others. Also describes their discovery
  method in more procedural detail (see "Method" below), and points to
  two resources: **TeloBase** (telomere motif database,
  http://cfb.ceitec.muni.cz/telobase/ - unreachable from this farm's
  network, proxy refused the connection when checked 2026-09-18; check
  manually before assuming a clade's motif is unknown) and **GERONIMO**
  (their own tool for automated CM-based ncRNA searches across many
  genomes with taxonomic/synteny context - not yet investigated whether
  it's usable here, worth following up before building redundant
  tooling).

## Why "beyond Hymenoptera/Lepidoptera" splits into two different problems

1. **Apply already-published TR models to DToL genomes nobody's searched
   yet** (vertebrates, fungi, etc. - biogenesis classes already
   characterized in the literature, e.g. via Podlevsky et al.'s
   Telomerase Database or Rfam). Engineering + reference-curation work,
   not discovery - out of scope for this `clade_discovery/` effort,
   which is specifically about (2).
2. **Genuine de novo discovery in clades with no TR characterization
   anywhere.** This is what `clade_discovery/` is for.

## Chosen target clades (decided 2026-09-18)

Picked by the user from a shortlist grounded in actual DToL taxonomic
depth (congeneric species count is the main feasibility lever - see
"Method" below for why):

| clade | DToL depth (directory-listing estimate) | **DToL depth (species with a FINISHED curated assembly - the number that actually matters)** | literature signal |
|---|---|---|---|
| **Araneae** (spiders) | `Clubiona` (7 spp), `Araneus` (5), several more 4-sp genera | **18 species, 22 genera matched**; best congeneric depth: `Araneus` (3), `Clubiona` (3), `Philodromus` (2) | Explicit telomerase-loss candidate, Fajkus & Fajkus 2025 citation [43] |
| **Myriapoda** (centipedes+millipedes) | `Lithobius` (9 spp - looked like the deepest non-Hymenoptera/Lepidoptera genus on the farm) | **Only 4 species total, 1 per genus - zero congeneric pairs currently available.** Most of `Lithobius`'s 9 directory entries are still mid-pipeline (no finished `primary.fa.gz` yet). This materially changes feasibility - see "Revised assessment" below. | Generic Ecdysozoa call-out, no specific prior claim |
| **Coleoptera** (beetles) | `Aphodius` (15), `Pterostichus`/`Longitarsus` (9 each) | **166 species, dozens of genera** - `Amara` (7), `Chrysolina`/`Cantharis` (6 each), `Pterostichus`/`Anaspis`/`Agonum` (5 each), many more with 3-4. **Turns out to be the deepest, most tractable of the three, despite having no specific literature hook.** | No specific literature hook found - speculative, included at user's request |
| **Nematoda** | 4 species (directory count = finished-assembly count here) | Same: 4 species, 2 in genus `Caenorhabditis`, 2 unrelated genera | Most explicitly named unknown in both papers, specific loss claim for *Meloidogyne* |

**Revised assessment (2026-09-18) after building the real assembly
lists**: directory *existence* under `/data/tol/data/darwin/` massively
overstates usable depth - `any_latest_assemblies.bash` only returns
species with a completed `primary.fa.gz`, and most DToL species
directories are still mid-curation-pipeline. Real depth ranking is now
**Coleoptera > Araneae >> Myriapoda ≈ Nematoda**. Myriapoda, despite
looking like the second-best pick from directory counts, currently has
no congeneric pairs at all in finished form - the reciprocal-filtering
step (Method step 4) has nothing to work with yet for this clade. Not
abandoning it, but it's now a "revisit once more assemblies finish
curating" clade rather than a pilot to build the discovery pipeline on
right now.

Nematoda is structurally different from the other three: DToL alone
doesn't have enough congeneric depth for the reciprocal-filtering trick
(see Method step 4) to work - would need genomes pulled in from NCBI,
a real extra step (data acquisition, storage, and the SUSPECT-screening
safety net built for DToL's specific curated-directory conventions
won't directly apply to NCBI assemblies). Treat as a later phase once
the DToL-native pipeline (Araneae/Myriapoda/Coleoptera) is built and
proven, not a first target.

## Method (generalized from Fajkus & Fajkus 2025's Figure 2 description)

1. **Genus/clade selection** - needs congeneric depth (multiple closely
   related, well-assembled genomes) AND, ideally, telomere-motif
   *variation* within the group. A clade with one static universal motif
   gives candidate-filtering nothing to work with - this is explicitly
   stated in the review: "focusing on species with telomere sequence
   variability can significantly reduce the number of potential
   candidates." We don't have per-species motif data for any of these
   clades yet - that's Stage 0 below, and real prioritization within
   these three clades (which genus/family to focus deepest effort on)
   should wait for it.
2. **Telomere motif discovery** per species - `tidk explore` (de novo,
   no assumed motif needed). This is Stage 0.
3. **Candidate TR-locus prediction** - template-shaped regions (any
   rotation of the discovered telomere repeat, ±1nt) extended with
   genomic context, in one reference genome per clade.
4. **Reciprocal blastn filtering** across the closely-related genus/
   family set - strips the (typically hundreds of) candidates down to
   the few that are shared/conserved. This step is *why* congeneric
   depth is the main feasibility lever - a lone, unrelated genome gives
   nothing to filter against.
5. **Cross-clade homology expansion** - blastn against progressively
   more distant relatives to find the survivors' taxonomic boundary.
6. **Covariance model (CM) construction** via Infernal - captures
   secondary structure, not just sequence; needed because TR sequence
   conservation is low even where structure is conserved. New territory
   for this project (the parent plant pipeline only ever used linear
   HMMs/nhmmer, not CMs/Infernal).
7. **Iterative refinement** - expand the CM outward per sub-clade.
8. **Structural validation** - MSA-guided secondary structure prediction
   (pseudoknot/PK, template boundary element/TBE, P1c stem-loop) to
   argue a candidate is a genuine TR, not incidental sequence.
9. *(Not planned here - wet-lab)*: RNA-seq (total RNA-seq with rRNA
   depletion specifically - the 2023 paper found poly-A-selected
   libraries produce false negatives for this TR class), northern blot,
   TRAP assay, FISH. Out of scope for this computational project; noting
   for completeness/honesty about the evidentiary ceiling here.

An important promoter-based alternative anchor (used by Fajkus for
early-diverging/more divergent clades, in addition to template
matching): scan for the type-3 snRNA promoter (USE + TATA box,
characterized via Infernal searches against Rfam's U1-U6/MRP/7SK CMs,
then MEME on the ~200nt extended regions) and require co-occurrence
with a template-like region within 200nt downstream. Useful when
telomere motif itself is uncertain/unconfirmed - not our primary route
for Araneae/Myriapoda/Coleoptera since we can discover motifs directly
via tidk, but worth keeping in reserve, and is essentially required if
we ever attempt a clade with no confirmed telomere motif in hand.

## What's built so far (this session)

- `src/classify_genera_ncbi.py` - classifies genus names via NCBI's
  Datasets API `group_name` field (e.g. "spiders", "centipedes",
  "beetles") so target-clade genus lists are built from real taxonomy,
  not memory/guesswork. Cached to `inputs/genus_ncbi_classification.tsv`.
- `src/build_clade_assembly_list.bash` - builds a clade-scoped
  species/path assembly list, reusing the parent pipeline's proven
  organelle/cobiont SUSPECT screening (naming-pattern + size-floor +
  plastid/mito keyword checks - see
  `../notes/organelle_and_cobiont_assembly_pitfall.md` for why this
  matters, a real DToL gotcha already hit twice in the plant pipeline).
- `src/run_tidk_clade.bash` - Stage 0 telomere motif discovery via
  `tidk explore`, LSF-submitted, per clade.

## Stage 0.5: positional confirmation (2026-09-21)

Stage 0's top genome-wide repeat is not necessarily telomeric (satellite/
centromeric repeats routinely outrank the true, much shorter telomeric
array by raw count) - so before touching candidate-locus prediction,
ran the actual positional check: `tidk search` restricted to chromosome
termini, for each species' top-3 Stage-0 candidates (not just #1).
Script: `src/confirm_repeat_positional_clade.py` (per-species worker,
deletes its own 10-20MB+ raw window dump immediately after analysis -
this farm's disk quota is a real, binding constraint, already hit once
this project), `src/run_confirm_positional_clade.bash` (LSF submission,
one job per species), `src/merge_positional_confirmations.bash`
(combines the per-species summaries - kept one-file-per-species during
the run specifically to avoid concurrent-write corruption across ~190
parallel jobs). Output: `outputs/positional_confirmations.tsv`.

**Result - Coleoptera is a clear, strong hit; Araneae/Myriapoda much
weaker:**

| clade | species checked | species with >=1 confirmed terminal motif |
|---|---|---|
| Coleoptera | 166 | **106 (64%)** |
| Araneae | 18 | 4 (22%) |
| Myriapoda | 4 | 1 (25%, but n=4) |

**Coleoptera's confirmed motifs converge cleanly on a 5bp-periodic
family** - the large majority are `AACCT` or exact multiples/rotations
of it (`AACCTAACCT`, `AACCTAACCTAACCT`, `ACCTGACCTG...` - a frame-shifted
read of the same unit), 106 species, no wet-lab input, no assumed motif
going in. This is consistent with the well-established ancestral insect
telomere repeat family (`TTAGG`-type) and is a genuinely independent,
dataset-wide computational rediscovery of it - a strong sanity check
that Stage 0 + Stage 0.5 work correctly end-to-end. **One real minority
variant stands out**: 5 species confirmed with `AACCCAACCCAACCC...`
(`AACCC`, not `AACCT` - a one-base swap in the repeat unit) - this is
exactly the kind of within-clade motif *variability* the whole
discovery method (see "Method" above, step 1) depends on to make
candidate-TR-locus prediction tractable. Worth identifying which 5
species these are and treating them as the seed for reciprocal-
filtering, mirroring exactly how Fajkus used `B. sylvestris`'s point
mutation against the other *Bombus* species.

**Araneae's 4 confirmed hits look qualitatively different and less
trustworthy** - `AAAAAAC`, `AATATATATATATATAT`, `AAAAAATT`, `ATCGC` -
none resemble each other, and several look more like generic AT-rich
low-complexity sequence than a genuine repeat family. Possible
explanations, not yet distinguished: (a) spider telomere motifs really
are this heterogeneous/poorly conserved (plausible - Araneae's
telomere biology is much less characterized than insects'), or (b) some
of these are false positives from unresolved/low-complexity sequence
clustering near contig ends in draft assemblies (assembly-artifact
enrichment, not real telomere biology - the same class of confound
`check_gene_synteny_gap.py`'s N-content check was built for in the
parent plant pipeline). **Needs an N-content-near-terminus sanity check
before trusting any of these 4** - not yet done.

## Step 3, take 2: promoter-based anchoring (2026-09-22)

Naive template-matching for candidate TR loci (first version of
`predict_tr_candidate_loci.py`, tested on the 3 confirmed-`AACCC`
*Aphodius* species) turned out not to scale down to the literature's
"needle in a haystack" tractability even after correcting a specificity
bug (fixed copy-count threshold vs. one that scales with unit length -
see the script's docstring): `Aphodius_granarius` alone produced
460,440 candidates even after the fix. Short (5bp) repeat units are
just fundamentally less specific than the literature's 11-13bp units -
more tandem copies would help precision but risks conflating candidates
with the real telomeric/interstitial array. Decided (user's call) to
build the promoter-based anchor instead of brute-forcing reciprocal
blastn on ~460k candidates - **explicitly framed as reusable machinery
for other clades, not a one-off fix for Coleoptera specifically.**

**Built:**
- `inputs/rfam/snRNA_promoter_models.cm` - Rfam's U1/U2/U3/U4/U5/U6/
  U6atac/RNase_MRP/Arthropod_7SK covariance models (the exact 9 families
  Fajkus's Hymenoptera paper used), fetched from the full Rfam.cm
  release via `cmfetch`, then `cmpress`-indexed. ~2MB total - the large
  intermediate download (44MB Rfam.cm) was deleted after extraction.
- `src/run_snrna_search_clade.bash` - LSF-submitted `cmsearch` against
  these 9 models, per species. ~15 min wall-clock (8 CPUs) for a
  ~300Mb beetle genome in testing.
- `src/extract_snrna_promoters.py` - pulls confident hits (e<1e-5, a
  stricter cutoff than Infernal's own default inclusion flag, which was
  too permissive for this purpose - 417 "included" U5 hits vs 36 at
  e<1e-5 for one species) and extracts the upstream flanking region
  (200bp default), handling genomic strand correctly (a '-strand hit's
  "upstream" is the higher-coordinate flank, reverse-complemented back
  to transcript-relative orientation).
- `src/run_mafft.bash` - this farm's top-level `/software/team301/mafft`
  wrapper is broken (bash errors from unset internal path variables,
  produces empty output); worked around by calling the underlying
  binary directly with `MAFFT_BINARIES` set explicitly. Worth knowing
  about beyond clade_discovery if the parent plant pipeline ever needs
  MAFFT again (README lists it as a dependency).
- `src/build_promoter_consensus.py` - per-column consensus + a
  conservation-fraction track from a MAFFT alignment of extracted
  promoter regions, plus a scan for the best-conserved short window
  (the TATA-box candidate).

**Result, pilot on `Aphodius_granarius` (the AACCC-confirmed species):**
142 confident (e<1e-5) snRNA hits, all 9 families represented (U1: 16,
U2: 36, U3: 14, U4: 4, U5: 36, U6: 32, U6atac: 1, RNase_MRP: 2,
Arthropod_7SK: 1). Promoter-region alignment + consensus found a
**clear, independently-replicated TATA-box-like element** (`TATATA`
core, best single window `GTATATAT`, ~69% mean conservation) at a
consistent position near the transcript-proximal end of the extracted
region - critically, this shows up independently across genuinely
different, unrelated gene families (U2, U3, U4, U5, U6, U6atac,
RNase_MRP, Arthropod_7SK all show the same core pattern individually,
not just as an artifact of one family's internal tandem-array
redundancy - checked explicitly by inspecting raw, ungapped sequence at
the conserved column across families, not just trusting the pooled
consensus number). This is a real, clean, unseeded rediscovery of a
canonical type-3 RNAPIII promoter TATA box, directly analogous to what
Fajkus describes doing by hand - found here with no assumed sequence.
No clear second (USE/PSEA) element identified yet in this quick pass -
worth a closer, more careful look, but the TATA signal alone is already
a much stronger anchor than bare template-matching.

**Not yet done:**
- Run the same promoter search on the other 2 *Aphodius* species (and
  ideally `Geotrupes_spiniger`) to check the TATA consensus is
  consistent across species, not a one-genome fluke.
- Use the characterized promoter signature (candidate TATA-box-like
  sequence within ~200bp of a genomic position) AS the new specificity
  filter for candidate TR loci, in place of/alongside bare template
  matching - i.e. redo step 3 requiring BOTH a promoter-like signature
  AND a nearby template match, which should be a much stronger joint
  filter than either alone and should get candidate counts back down to
  something reciprocal-blastn can actually handle.
- Identify a candidate USE/PSEA element more carefully (may need a
  cleaner, less redundant sequence set - collapsing near-identical
  same-locus tandem copies before alignment would reduce bias from
  within-array redundancy).

## Not yet built (needs real Stage 0 data first to design sensibly)

- Candidate TR-locus prediction (step 3) - parameters (window size,
  rotation-matching) should be informed by what the motifs actually
  look like per clade once Stage 0 completes.
- Reciprocal blastn filtering (step 4).
- Infernal/CM construction (step 6) - new tooling dependency, not yet
  confirmed available on this farm (`which infernal`/`cmbuild`/`cmsearch`
  not yet checked).
- Structural validation (step 8).

## Open items to resolve before going further

- Check TeloBase manually (network-blocked from this environment) for
  existing motif records in Araneae/Myriapoda/Coleoptera - could save
  real Stage 0 effort for species already characterized.
- **RESOLVED (2026-09-18): Infernal is already installed** at
  `/software/team301/infernal-1.1.4-linux-intel-gcc/binaries/`
  (`cmbuild`, `cmsearch`, `cmalign`) - step 6 is reachable with no new
  software install.
- **RESOLVED (2026-09-18): GERONIMO is real, public, and directly
  relevant.** Open source (github.com/amkilar/GERONIMO, Fajkus group),
  published in GigaScience (10.1093/gigascience/giad080), Snakemake +
  Infernal based, purpose-built for exactly steps 4-7 (CM-based
  homology search across many genomes with taxonomic context, Excel/
  heatmap summaries). Tested at "hundreds of genomes" scale - matches
  our need. Needs evaluation before committing to build our own
  steps 4-7 from scratch: (a) it auto-downloads genomes from NCBI by
  query by default - need to check whether/how it accepts local genome
  files instead, since we specifically want DToL's curated assemblies,
  not whatever happens to be on NCBI for these UK species; (b) needs
  either an Rfam CM or a user `.stk` alignment as the seed model - we'd
  supply our own template-matched candidate alignment as the seed for
  a clade with no existing Rfam TR model. Strong candidate to adopt
  rather than reimplement - next concrete step once Stage 0 data is in
  hand.
