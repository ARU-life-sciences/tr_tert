"""Score conserved template-bearing candidates with TR models.
Usage: rank_candidates.py <fasta> <candidate ids file> <cm> <hmm library> <out prefix>
Candidates are seq:a-b:strand (covariation.py); each is extracted +/- 200 nt
in TR orientation and searched with cmsearch and nhmmer (E <= 10)."""
import sys, gzip, subprocess
sys.path.insert(0, '/lustre/scratch122/tol/teams/blaxter/users/mb39/tol_telomeres/shifts/tr_copies/src')
from copies import read_genome, rc
fa, ids, cm, hmm, out = sys.argv[1:6]
want = [l.strip() for l in open(ids) if l.strip()]
G = read_genome(fa)
with open(out + '.fa', 'w') as fh:
    for c in want:
        s, r, st = c.rsplit(':', 2); a, b = (int(x) for x in r.split('-'))
        seq = G[s][max(0, a - 201):b + 200]
        fh.write(f'>{c}\n{rc(seq) if st == "-" else seq}\n')
subprocess.run(['/software/team301/infernal-1.1.4-linux-intel-gcc/binaries/cmsearch', '--cpu', '4', '-E', '10',
                '--tblout', out + '.cm.tbl', cm, out + '.fa'], check=True, stdout=subprocess.DEVNULL)
subprocess.run(['/software/team301/hmmer-3.4/src/nhmmer', '--dna', '--cpu', '4', '-E', '10',
                '--tblout', out + '.hmm.tbl', hmm, out + '.fa'], check=True, stdout=subprocess.DEVNULL)
best = {}
for f, ecol, kind in ((out + '.cm.tbl', 15, 'cm'), (out + '.hmm.tbl', 12, 'hmm')):
    for l in open(f):
        if l.startswith('#'): continue
        x = l.split(); e = float(x[ecol])
        if x[0] not in best or e < best[x[0]][0]: best[x[0]] = (e, kind, x[2])
for c, (e, kind, model) in sorted(best.items(), key=lambda x: x[1][0])[:15]:
    print(c, f'{e:.1e}', kind, model, sep='\t')
print(f'{len(best)}/{len(want)} candidates with any model hit at E<=10', file=sys.stderr)
