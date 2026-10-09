# TR templates against the ToL telomere calls (2026-10-09)

Known TR models searched across whole clades, and each species' TR
template tested against its tidk telomere call (`tr_hym/src/template_check.py`:
a template supports a unit if it exactly matches >= len(unit) + 2 nt of
its tandem). Shared scripts in `tr_hym/src/`, run per clade with
`CLADE_DIR` and `MODEL_NAME` (and `SEARCH=cmsearch` for an Infernal CM).

| clade | model | species | TR found | strong call: template agrees |
|---|---|---|---|---|
| Hymenoptera (`tr_hym`) | 14 family HMMs from 256 published TRs (Fajkus et al. 2023, Table S1) | 298 | 182 | **110/115** |
| Lepidoptera + Trichoptera (`tr_lep`) | 20 family HMMs from 708 TRs (Fajkus et al. 2023, Table S4) | 911 | 910 | **800/803** |
| vertebrates (`tr_vert`) | Rfam RF00024 (Telomerase-vert CM) | 286 | 282 | **177/178** |
| **total** | | 1,495 | | **1,087/1,096 (99%)** |

Weakly called species agree too: Lepidoptera 107/107, vertebrates 97/101.

**The exceptions are mostly tidk miscalls, not template disagreements.**
In *Sabra harpagula* (tidk ACCTC), *Mechanitis lysimnia* (AAACC) and
*Aphomia sociella* (AACC) the template reads AACCT, the Lepidoptera
repeat; in *Ficedula hypoleuca* tidk called a 16-nt satellite and the TR
does not support it. The others: *Agrypon* (tidk called a compound of its
template unit), three *Andrena* (a genus with 7 units in 8 species), and
*Alloplasta piceator* (no template-like match). So the TR results are also
an independent validation of the telomere calls.

No TR hit: Hymenoptera 74 strongly called species, mostly sawflies and
Ichneumonidae (few or no published TRs); Lepidoptera only *Adscita geryon*
(Zygaenidae).

Data problems found on the way:
- The curated assemblies of *Sebastes mentella* (fSebMen2.1) and *Myodes
  glareolus* (mMyoGla1.1) contain runs of NUL bytes inside sequence lines;
  `tr_hym/src/worker_search_v2.bash` strips them.
- nhmmer needs `--dna`: T2T assemblies can start with kilobases of G-free
  telomere, which defeats its alphabet guess.

De novo: `denovo/README.md` (template co-variation; works within *Bombus*,
not yet elsewhere).
