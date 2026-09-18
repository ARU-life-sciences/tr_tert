# Epilobium hirsutum: TERT is NOT genuinely absent - it's present on hap1, missing specifically from hap2

**Status: REVISED 2026-09-18 - original "candidate TERT loss" finding
does not hold at the species level.** Checking the alternative haplotype
(prompted by this note's own previously-listed caveat, "only one
haplotype checked") found a strong, unambiguous TERT hit on
*E. hirsutum* haplotype 1, at the exact syntenic locus where haplotype 2
(the assembly used throughout the rest of this pipeline, and in the
original version of this finding below) shows nothing. This is now a
candidate **haplotype-specific** absence - real and still worth
reporting, but a materially different, more modest claim than "TERT
genuinely lost from this species." See "Revision" below for the new
evidence; the original hap2-only analysis is kept beneath it for the
record.

## Summary

*Epilobium hirsutum* (great willowherb, Onagraceae) shows no confident
TERT (telomerase reverse transcriptase) hit by any of 3 HMM profiles or
by tblastn against **haplotype 2** of its genome (`drEpiHirs1.hap2.1`,
the assembly this pipeline tracks and searches by default) - the same
surface pattern as `Viscum_album` (headline #1). But checking
**haplotype 1** (`drEpiHirs1.hap1.1`) directly finds a strong,
multi-exon TERT hit (all 3 HMM profiles, e-values to 4.9e-48; tblastn,
334 hits, e-values to 1.5e-111) at `SUPER_13:~15,367,500-15,369,250`,
confirmed by two independent synteny comparisons to be the same physical
locus as hap2's gap (see "Revision" below). So TERT is not missing from
this species/individual overall - it's missing from one haplotype
assembly and present on the other.

Its TR (telomerase RNA) gene is, separately, clearly present (3 loci,
all sharing an identical core template, `AAACCCTAAACCT`) and was already
independently confirmed elsewhere in this investigation (the
TR-identical control group built for headline #6) to genuinely match
real telomeric sequence in an exhaustive, base-pair-resolved,
telomere-proper-restricted walk - sparse (3 confirmed occurrences
dataset-wide) but real and non-zero. That part of the original finding
is unaffected by this revision.

## Background

The standard pipeline searches every DToL species for TERT via 3 HMM
profiles (`catalytic_subunit`, `isoformx1`, `tert1` - `run_tert_hmms.bash`)
followed by tblastn confirmation against curated plant TERT references
(`run_tblastn_tert_bsub.bash`). `Epilobium_hirsutum` returns zero hits
from all four searches:

```
outputs/tert_tbls/Epilobium_hirsutum/*.tbl   - all 3 HMM profiles: 0 hits
outputs/tert_tblastn/Epilobium_hirsutum.tbl  - 0 bytes
```

Its TR gene search (`outputs/tr_tbls/Epilobium_hirsutum.tbl`), by
contrast, is a normal, non-empty, multi-locus hit: 3 loci, all sharing
one identical 13bp core template (`AAACCCTAAACCT`) -
`core_template_sets.tsv` row `Epilobium_hirsutum 13 3 1 AAACCCTAAACCT`.
This made it eligible for, and it was included in, the TR-identical
control group built for headline #6
(`outputs/tr_repeat_correlation/identical_set_array_structure.tsv`),
where the same exhaustive, telomere-proper-restricted walk confirmed the
TR core genuinely present in the real telomere (`proper_dom_n=3` - a low
but non-zero, real confirmed count; the candidate second/minor variant
tested alongside it in that analysis was not confirmed present, unrelated
to the point here). So the TR side of telomerase is not in question;
only TERT.

The obvious alternative explanation for "no hit" is assembly
incompleteness (the same caveat already flagged, unresolved, for
Viscum), not genuine absence - draft genome assemblies can have real
gaps, and a single-copy gene sitting in one is indistinguishable from a
genuinely deleted one by a simple hit/no-hit search alone. This needed a
positive, independent check, not just "we didn't find it."

## Evidence: synteny check against a congener with a confident TERT hit

`Epilobium_ciliatum` (same genus) has an unambiguous, strong TERT hit -
consistent multi-exon matches across all 3 HMM profiles and tblastn, at
e-values down to 1e-36, all localised to one region of `SUPER_16`
(~2,602,800-2,609,900). That rules out "our TERT reference library
doesn't cover this lineage" - a close relative is trivially found.

`src/check_gene_synteny_gap.py` extracts a 20kb window around this
region in `E. ciliatum` (`SUPER_16:2,595,000-2,615,000`, generously
padded on both sides of the gene) and blastn's it against `E. hirsutum`'s
full genome (e-value <= 1e-10):

```bash
python3 src/check_gene_synteny_gap.py Epilobium_ciliatum SUPER_16 2595000 2615000 \
  Epilobium_hirsutum outputs/synteny_checks
```

Result:

- **A confident, single-locus anchor exists**: query positions
  3,794-6,296 (i.e. `E. ciliatum` genome coordinates ~2,598,800-2,601,300,
  just upstream of the TERT gene) match `E. hirsutum SUPER_16:3,043,979-
  3,047,576` at 87-93% identity, colinear, best e-value 0.0 - a real,
  single-locus syntenic anchor, not scattered repetitive-element noise
  (repetitive-element hits, by contrast, scatter across a dozen-plus
  different scaffolds at once - seen clearly flanking this anchor on
  both sides in the raw hit table, and correctly not treated as
  informative synteny).
- **The TERT gene body itself has zero hits**: query positions
  6,883-15,353 (~8.5kb, comfortably spanning the entire TERT ORF region
  based on `E. ciliatum`'s own coordinates) returns **no blastn hit
  anywhere in `E. hirsutum`'s genome**, at any identity, even at the
  permissive 1e-10 threshold used for the whole search.
- **Not an assembly gap at this locus**: the anchor sits on
  `SUPER_16`, a 10,471,454bp scaffold, at position ~3.04Mb - 29% of the
  way through, nowhere near either end. A 30kb window centred on the
  anchor (`SUPER_16:3,030,076-3,060,076`) contains **0% N's** - fully
  resolved, real, sequenced DNA, not unresolved/gapped assembly.

So the genomic neighbourhood is present, confidently localised to one
specific place, and completely resolved (in hap2) - but whatever
sequence actually occupies the interval where TERT should be, it isn't
detectably TERT there, at a locus where the same search finds TERT
easily in a close relative sharing 87-93% flanking identity. **This is
where the original version of this note stopped - see "Revision" below
for what checking hap1 changed.**

## Revision (2026-09-18): checking haplotype 1 finds TERT is present there

The Caveats section below (before this revision) already flagged the
obvious next check as unfinished: only hap2 had been searched, and the
congener comparison used a different haplotype (`E. ciliatum` hap1)
than `E. hirsutum` (hap2) - "worth confirming `E. hirsutum`'s hap1 (if
available) agrees before treating this as settled." It doesn't agree.

**TERT is present on `E. hirsutum` hap1** (`drEpiHirs1.hap1.1.primary`,
run directly rather than through the pipeline's tracked path, which
points at hap2):

- All 3 HMM profiles hit the same locus, `SUPER_13:~15,367,500-15,370,550`,
  as a clean multi-exon/multi-segment match (consistent with a real,
  intron-interrupted gene in genomic DNA, the same pattern seen for
  `E. ciliatum`'s confirmed hit): `isoform_X1` e-value 4.9e-48 (covering
  hmm columns 221-1690 of 1813, ~81% of the model), `catalytic_subunit`
  e-value 7.7e-25, `tert1` e-value 9.9e-24.
- **tblastn against the same locus: 334 hits**, e-values to 1.5e-111,
  ~46% identity, against multiple plant TERT references (`KAI9397291.1`,
  `XP_024452485.2`, and others) - far stronger than the borderline
  e-value range (down to 1e-4) the standard pipeline already treats as a
  confident call.
- **0% N's** in the region spanning the hit - not a hap1 assembly gap
  either.
- Position: 87% of the way through a 17.6Mb scaffold - not
  scaffold-terminal.

**This hap1 locus is the same physical locus as hap2's gap**, confirmed
two independent ways:
1. The original congener-based anchor (`E. ciliatum` vs `E. hirsutum`
   hap2) landed at `SUPER_16:3,043,979-3,047,576`.
2. Blasting a window around the new hap1 TERT hit directly against hap2
   (bypassing `E. ciliatum` entirely) finds a small cluster of
   high-identity (87-93%), non-repetitive hits at `SUPER_16:3,043,979-
   3,046,585` - **the start coordinate matches the first anchor to the
   exact base pair** (3,043,979), and sits ~2.3-3.9kb downstream of the
   hap1 gene body, i.e. clearly the same flanking sequence, reached via
   two unrelated routes. (The rest of the ~38kb window around the hap1
   gene is highly repetitive - scattered 75-93%-identity hits across 87
   different hap2 scaffolds - and was correctly not treated as
   informative synteny, same logic as the original check's repetitive-
   noise exclusion.)

So: real, multi-method, high-confidence TERT gene on hap1; a confirmed
absence on hap2; and independent proof both loci are the same place in
the genome, not two unrelated regions. The original "genuine TERT loss"
framing does not hold at the species/individual level - the gene is
there, just not on the haplotype this pipeline happens to track.

## Interpretation

**Revised.** This is now a candidate **hap2-specific TERT absence**
within one individual, not a species-level loss. Two explanations
remain open, and this analysis can't distinguish them:

1. **Genuine hemizygous/allelic deletion** - this individual's hap2
   chromosome copy really has lost the TERT gene body (a real structural
   difference between its two haplotypes), while hap1 retains an intact
   copy. Biologically plausible and would itself be a notable finding
   (a real, sequenced case of loss-in-progress at one allele), but
   distinct from and weaker than "TERT lost from the species."
2. **A hap2-specific assembly artifact** localised just to this gene's
   span - despite passing every check this analysis can run (flanking
   region intact, confident single-locus anchor, 0% N's, safely
   mid-scaffold), a small, gene-body-sized collapse or deletion
   specifically in hap2's assembly (as opposed to broad incompleteness,
   which was already ruled out) isn't something a synteny/N-content
   check alone can fully exclude - it would need read-level coverage
   evidence across the hap2 locus, which wasn't available for this
   check.

Either way, the TR side of the story is unaffected: TR is present,
multi-copy, and independently confirmed to template the real telomere.
What's no longer supported is the "both components effectively gone in
this species" framing this note originally proposed as a parallel to
Viscum - that comparison should be dropped.

## Caveats

- **Whether hap2's absence is real biology or a hap2-specific assembly
  artifact is unresolved** (see Interpretation above) - this is now the
  central open question, not "is TERT missing from the species," which
  is resolved (no, it isn't).
- **Single congener comparison, one direction** (for the original
  hap2-vs-`E. ciliatum` check) - doesn't rule out a large-scale
  rearrangement having moved rather than deleted the gene, though this
  is now largely moot since hap1 already shows the gene at the
  expected, syntenically-confirmed locus.
- **No wet-lab validation** - this remains purely comparative-genomic
  evidence.
- **The exact gene-body boundaries** are still defined relative to
  reference coordinates (`E. ciliatum`'s gene model, or the hap1 hit
  coordinates), not an independently annotated `E. hirsutum` gene model.
- **No literature check** has been done for this species or genus
  (Onagraceae), including specifically for prior reports of
  intraspecific/allelic telomerase-component variation.
- **The standard pipeline still tracks `E. hirsutum` hap2 as this
  species' assembly** (`inputs/dtol_plant_paths.txt`), so
  `tert_absence_summary.py`'s dataset-wide completeness check (and any
  other pipeline output) will keep reporting this species as
  TERT-absent unless/until that's revisited - a live inconsistency to
  be aware of, not yet resolved here. Worth checking, at some point,
  whether other species that surface as candidate absences have the
  same haplotype-specific pattern before trusting a hap2-only "missing"
  call as species-level.

## Reproducing

```bash
# Confirm the TR side is solid (already established elsewhere in this repo):
grep -c "." outputs/tr_tbls/Epilobium_hirsutum.tbl   # non-empty TR hit
grep "^Epilobium_hirsutum" outputs/tr_repeat_correlation/core_template_sets.tsv
grep "^Epilobium_hirsutum" outputs/tr_repeat_correlation/identical_set_array_structure.tsv

# Confirm TERT absence on hap2, across all 4 standard searches
# (this is what the pipeline tracks/searches by default):
for f in outputs/tert_tbls/Epilobium_hirsutum/*.tbl; do grep -vc "^#" "$f"; done   # all 0
wc -l outputs/tert_tblastn/Epilobium_hirsutum.tbl                                  # 0 bytes

# Confirm the congener has a real hit (the positive control):
head -5 outputs/tert_tblastn/Epilobium_ciliatum.tbl

# Original synteny check (E. ciliatum vs E. hirsutum hap2):
python3 src/check_gene_synteny_gap.py Epilobium_ciliatum SUPER_16 2595000 2615000 \
  Epilobium_hirsutum outputs/synteny_checks

# Revision: TERT HMM + tblastn run directly against E. hirsutum hap1
# (not the pipeline's tracked hap2 path) - see run_tert_hmms.bash /
# worker_tblastn_tert.bash for the exact commands these mirror:
HAP1=/data/tol/data/darwin/dicots/Epilobium_hirsutum/assembly/curated/drEpiHirs1.1/drEpiHirs1.hap1.1.primary.fa.gz
gunzip -c "$HAP1" > /tmp/Epilobium_hirsutum_hap1.fa
for hmm in isoform_X1_tert catalytic_subunit_tert tert1; do
  /software/team301/hmmer-3.4/src/nhmmer --dna --cpu 8 \
    --tblout "hap1_${hmm}.tbl" "inputs/clusters/${hmm}.hmm" /tmp/Epilobium_hirsutum_hap1.fa
done
bash src/worker_tblastn_tert.bash Epilobium_hirsutum_hap1 "$HAP1"

# Confirm the hap1 hit and the hap2 gap are the same physical locus,
# independent of the E. ciliatum comparison above - extract a window
# around the hap1 hit and blast it directly against hap2:
esl-sfetch -c "15352000..15390000" /tmp/Epilobium_hirsutum_hap1.fa SUPER_13 > hap1_window.fa
gunzip -c /data/tol/data/darwin/dicots/Epilobium_hirsutum/assembly/curated/drEpiHirs1.1/drEpiHirs1.hap2.1.primary.fa.gz > /tmp/Epilobium_hirsutum_hap2.fa
makeblastdb -in /tmp/Epilobium_hirsutum_hap2.fa -dbtype nucl -out hap2_db
blastn -query hap1_window.fa -db hap2_db -evalue 1e-10 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  | awk -F'\t' '$2=="SUPER_16"'   # look for the ~3,043,979-3,046,585 cluster
```
