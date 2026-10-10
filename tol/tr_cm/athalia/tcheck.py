import gzip,sys
sys.path.insert(0,'/lustre/scratch122/tol/teams/blaxter/users/mb39/tol_telomeres/shifts/tr_hym/src')
from template_check import longest_tandem_match as m
COMP=str.maketrans('ACGTN','TGCAN')
for sp,f,hits in [('Athalia_cordata',sys.argv[1],[('SUPER_1',40834659,40834515)]),('Athalia_circularis',sys.argv[2],[('SUPER_2',23309504,23309364),('SUPER_4',5800804,5800681),('SUPER_1',28333578,28333708)])]:
    want={h[0] for h in hits}; seqs={}; name=None; buf=[]
    for line in gzip.open(f,'rt'):
        if line.startswith('>'):
            if name in want: seqs[name]=''.join(buf)
            name=line[1:].split()[0]; buf=[]
        elif name in want: buf.append(line.strip().upper())
    if name in want: seqs[name]=''.join(buf)
    for c,a,b in hits:
        lo,hi=min(a,b),max(a,b); s=seqs[c][lo-301:hi+300]
        if a>b: s=s.translate(COMP)[::-1]
        print(sp,c,a,b,'AACCT match',m(s,'AACCT'),'| CCTAACCTAAC in hit:', 'CCTAACC' in s)
