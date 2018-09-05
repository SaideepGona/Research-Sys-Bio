#!/usr/bin/env python
import pdb,sys,os

# Filters a bed file based upon peak size

f=open(sys.argv[1])
lf=f.readlines()
f.close()
lf=[item.strip().split() for item in lf]
out=[]

cut=int(sys.argv[2])
# pdb.set_trace()
for i in lf[1:]:
    try:
        istart=int(i[1])
        iend=int(i[2])
        # pdb.set_trace()
        span= max(iend,istart)-min(iend,istart)
        if span >cut:
            out.append(i)
    except:
        pass
# pdb.set_trace()
out=["\t".join(item) for item in out]
out="\n".join(out)

f=open(sys.argv[1]+"_filter.txt",'w')
f.write(out)
f.close()
