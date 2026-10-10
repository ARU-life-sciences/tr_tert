import gzip,subprocess,sys
H='/tmp_relic'
checks=[('Scambus_brevicornis','SUPER_6_HAP1',1335278,1336963),('Dolichomitus_mesocentrus','SUPER_10',20776216,20778266),('Acrodactyla_quadrisculpta','SUPER_5_HAP1',56150115,56187018)]
import csv
fa={r['species']:r['fasta'] for r in csv.DictReader(open('tr_hym/species.tsv'),delimiter='\t')}
for sp,seq,lo,hi in checks:
    s=[];on=False
    for line in gzip.open(fa[sp],'rt'):
        if line.startswith('>'):
            if on: break
            on=line[1:].split()[0]==seq; continue
        if on: s.append(line.strip())
    s=''.join(s); lo2,hi2=max(0,lo-5000),hi+5000
    open('tr_pimpl/relic.fa','w').write(f'>{seq}_{lo2}\n{s[lo2:hi2]}\n')
    out=subprocess.run(['/software/team301/hmmer-3.4/src/nhmmer','--dna','-E','10','--tblout','/dev/stdout','tr_pimpl/models/pimplinae_tr.hmm','tr_pimpl/relic.fa'],capture_output=True,text=True).stdout
    hits=[l.split() for l in out.splitlines() if l and not l.startswith('#') and len(l.split())>12 and l.split()[0].startswith(seq)]
    print(sp, f'{seq}:{lo2}-{hi2} ({(hi2-lo2)/1000:.1f} kb):', [(int(h[6])+lo2, int(h[7])+lo2, h[12]) for h in hits] or 'no hit at E<=10')
