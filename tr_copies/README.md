# What are the extra TR-like copies? (2026-10-09, preliminary)

`src/copies.py`: for every TR-like locus (E <= 1e-5, one per 2 kb) in 9
pimplines, the 30 Hymenoptera with most copies, and 20 single-copy
Lepidoptera (baseline), with 20 random positions per genome as background:
a poly-A run after the 3' end, a target-site duplication (direct repeat >=
10 nt either side), a TATA-like box 20-40 nt upstream, and a template.

| set | copies | poly-A | TSD >= 10 nt | TATA | template |
|---|---|---|---|---|---|
| Lepidoptera TRs | 23 | 0% | 0% | **74%** (random 22%) | 100% |
| pimplines | 38 | 13% (random 3%; 5/30 decayed copies) | 0% (random 9%) | 13% (random 17%) | 21% |
| multi-copy Hymenoptera | 438 | 0.2% (random 0.8%) | 19% (random 30%) | 30% (random 13%) | 54% |

- The TATA test works on full-length single-copy TRs (Lepidoptera 74% vs
  22%); pimpline hits cover only the 3' core, so their "upstream" window is
  not the promoter.
- No retrocopy signature: poly-A and TSDs at background (the TSD test is
  noisy in AT-rich genomes).
- `src/flanks.py` (copies' 400-nt flanks blasted against each other;
  clustering):

| species | copies | sequences | copies sharing flank with another copy |
|---|---|---|---|
| *Lasioglossum calceatum* | 52 | 22 | 2/52 |
| *Seladonia tumulorum* | 27 | 13 | 16/27 |
| *Andrena nitida* | 23 | 4 (12 on one) | 18/23 |
| *Perithous albicinctus* | 7 | 5 | 0/7 |
| *Paraperithous gnathaulax* | 6 | 5 | 0/6 |

Two patterns: dispersed TR-sized insertions with unique flanks
(pimplines, *Lasioglossum*), and copies inside a larger repeated unit
(*Andrena*, *Seladonia*: a TR-carrying repeat or segmental duplication).
Retroposition is not supported so far; the mechanism is open. Next: full
TR-length copy alignments (5' ends), TE annotation of the flanks, copy
ages (divergence between copies).
