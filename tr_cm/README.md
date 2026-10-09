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
