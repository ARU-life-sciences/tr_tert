"""Do a genome's TR-like copies share flanking sequence (parts of one
larger repeat, e.g. a TR-derived SINE), and are they clustered?
Usage: flanks.py <species> <fasta> <hits.fa>"""
import sys, subprocess, collections
sys.path.insert(0, '/lustre/scratch122/tol/teams/blaxter/users/mb39/tol_telomeres/shifts/tr_copies/src')
from copies import read_genome, rc
sp, fasta, hits = sys.argv[1:4]
loci = {}
for h in (l[1:].strip() for l in open(hits) if l.startswith('>')):
    loc, model, e = h.split('|'); s, r = loc.rsplit(':', 1); x, y = (int(v) for v in r.split('-'))
    if float(e) > 1e-5: continue
    k = (s, min(x, y) // 2000)
    if k not in loci or float(e) < loci[k][3]: loci[k] = (s, min(x, y), max(x, y), float(e), '+' if x < y else '-')
G = read_genome(fasta)
per = collections.Counter(v[0] for v in loci.values())
pos = sorted((v[0], v[1]) for v in loci.values())
gaps = [b[1]-a[1] for a, b in zip(pos, pos[1:]) if a[0] == b[0]]
print(f"{sp}: {len(loci)} copies on {len(per)} sequences; max per sequence {max(per.values())}; "
      f"spacings < 10 kb: {sum(g < 10000 for g in gaps)}/{len(gaps)}")
for side in ('left', 'right'):
    with open(f'/tmp/{sp}_{side}.fa', 'w') as fh:
        for i, (s, a, b, e, st) in enumerate(loci.values()):
            reg = G[s][max(0, a-500):b+500]
            if st == '-': reg = rc(reg)
            flank = reg[:400] if side == 'left' else reg[-400:]
            fh.write(f'>c{i}\n{flank}\n')
    out = subprocess.run(['/software/team301/ncbi-blast-2.16.0+/bin/blastn', '-query', f'/tmp/{sp}_{side}.fa',
                          '-subject', f'/tmp/{sp}_{side}.fa', '-evalue', '1e-10', '-outfmt', '6 qseqid sseqid length pident'],
                         capture_output=True, text=True).stdout
    pairs = {tuple(sorted(l.split()[:2])) for l in out.splitlines() if l.split()[0] != l.split()[1]}
    with_match = {x for p in pairs for x in p}
    print(f"  {side} 400-nt flanks: {len(with_match)}/{len(loci)} copies share flank sequence with another copy ({len(pairs)} pairs)")
