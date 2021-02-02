#!/usr/bin/python

from glob import glob
from sys import argv
from tools.base import range2list
from tools.bctab import *
from tools.fraMotif import *

'''This tool is used to predict the chemical shifts of  RNA imino group within A-form helix. Only support base pairs [GC,CG,UA,AU,GU,UG] for now.''' 

#===============================================================================
NUS = {'G':('N1','H1'),'U':('N3','H3')}

#===============================================================================


if len(argv) <2:
    print '\n\n\tUsage: genimino.py <sequence_bracket-dot_file> ',\
    '<residue_range[optional]>\n',\
    "\te.g. genimino.py tP5abc.seq 130-145,150-153,158-193\n",\
    "\tor genimino.py tP5abc.seq\n\n"
    exit(1)

fnseq = argv[1]
lines = open(fnseq).readlines()
seq = lines[0].strip()
bdstr = lines[1].strip()
CS = BCTab('./tools/NH.cs').NH
resis = []
if len(argv)==3:
    resis = range2list(argv[2])
    if len(resis)!= len(seq):
        print "\tPlease check the consistency betwwen residue range and sequence!"
        exit(1)
    
outl =[] 
x= []
y = []
peaks = []
for i,s in enumerate(seq):
    if s not in NUS:
        continue
    if resis:
        resi = resis[i]
    else:
        resi = i+1
    idx,mode = getMotif('triplet',i,seq,bdstr)
    if mode:
        nuN = NUS[s][0]
        nuH = NUS[s][1]
        bd = '-'.join(mode[:3])
        try:
            N = CS[bd][0]
            H = CS[bd][1]
        except KeyError:
            print '\n\n\tWarning: Only base pairs [GC,CG,UA,AU,GU,UG] can be predicted for now.\n',\
            '\tPlease check the bracket-dot line! Residue: %d'%resi
            continue
        outl.append('%1s\t%5d\t%5s\t%8.2f\n'%(s,resi,nuN,N))
        outl.append('%1s\t%5d\t%5s\t%8.2f\n'%(s,resi,nuH,H))
        x.append(H)
        y.append(N)
        peaks.append(s+str(resi))
fout = open('imino.tab','wt')
#fout.write('ResName Residue Nuclei Chemical shift\n')
fout.writelines(outl)
fout.close()
print "\n\n\tOutput: predicted chemical shifts have been saved in <imino.tab>!\n\n"
#pyl.figure(1,figsize=(14,9))
#NHSpec(x,y,111,xlimt=[9.6,16.2],ylimt=[141,167],peaks=peaks)
#pyl.savefig('NHspec.png')
#pyl.show()
