# Pimpline TRs (2026-10-09)

A Pimplinae TR profile HMM (`models/pimplinae_tr.hmm`, from the 7
validated pimpline TRs: *Ephialtes*, *Scambus*, *Dolichomitus*,
*Acrodactyla*, *Allomacrus*, *Pimpla*, *Apechthis*; 79-92% identity
between genera) searched in the 9 AAAGAACCT/AAAGAACCC ichneumonids with
no TR after model walking (whose Ichneumonidae model mixed divergent
subfamilies).

Every one of the 9 genomes has a strong TR homolog (E 1e-13 to 1e-18), but
also **3-7 TR-like copies on different chromosomes, usually only one with
a template** for its telomere. The best-scoring copy is often a decayed
one, so the template check picks the best copy with a supporting template
(`PICK_TEMPLATE=1 template_check.py`):

| result | species |
|---|---|
| validated (template, E <= 1e-10) | *Schizopyga*, *Zaglyptus* (AAAGAACCC), *Paraperithous*, *Dreisbachia*, *Cylloceria* |
| template-bearing copy, weaker hit | *Clistopyga* (E 1e-5), *Perithous divinator* (1e-8), *Collyria* (2e-8) |
| no template (best 10 nt) | *Perithous albicinctus* |

`relic.py`: the intervals in *Scambus*, *Dolichomitus* and *Acrodactyla*
syntenic with the *Ephialtes* TR hold no TR-like sequence even at E <= 10
(12-47 kb): the *Ephialtes* TR is an insertion at a new site, not an
ancestral copy decayed in the others.
