# Covariance models for TR (2026-10-09)

`src/build_cm.py` (MAFFT L-INS-i, RNAalifold consensus structure, cmbuild).

- Sawfly (Tenthredinidae, 28 DToL TRs): 157 columns, 26 consensus base
  pairs (two hairpins and an enclosing stem). The 9 published
  Tenthredinoidea TRs, aligned on their own, give the same fold.
- Ichneumonidae: no alignable consensus over all 69 TRs (0 base pairs);
  the 9 published true ichneumonid TRs give 8 base pairs.
- Structure does not carry across sawfly lineages: under the
  Tenthredinidae-trained CM, Tenthredinidae TRs score a median 85 bits,
  but *Athalia rosae* 2, *Cephus cinctus* 3, *Orussus* -12, and true
  ichneumonids -20 to +1. A family-trained CM does not reach the next
  family.
- `src/run_sawfly_cm.bash`: CM from all 39 Tenthredinoidea TRs (cmalign to
  the seed, cmbuild, cmcalibrate), searched in the 9 sawflies still without
  a TR -> `outputs/sawfly_cm/`.

### Result: the Tenthredinoidea CM on the 9 remaining sawflies

| species | family | CM hit | template (tidk unit) |
|---|---|---|---|
| *Allantus rufocinctus* | Tenthredinidae | E 6.9e-11 | AACCT 11 nt: validated |
| *Athalia circularis* | Athaliidae | E 1.7e-5 | AACCT 12 nt |
| *Athalia cordata* | Athaliidae | E 3.8e-4 | AACCT 12 nt |
| *Abia aenea* | Cimbicidae | E 2.1e-5 | AACCT 14 nt |
| *Arge ustulata* | Argidae | E 9.9e-5 | no AACCC template: not a TR |
| *Aproceros*, *Arge ochropus*, *Urocerus*, *Xiphydria* | | none | |

The *Athalia* CM hits are the same loci as the *A. rosae* profile hits
(SUPER_2:23.31 Mb, SUPER_1:40.83 Mb): two independent models plus a
template make them good candidates despite weak E-values. Still no TR in
Argidae (whose telomere is AACCC, a shift), Siricidae or Xiphydriidae.
