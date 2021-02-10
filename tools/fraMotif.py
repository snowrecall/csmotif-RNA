from glob import glob
from base import range2list

'''
Fragment RNA Motif according to secondary structure
Supported motif types: triplet , penta, base pair
'''
tp = {'19':['GC','CG'],'20':['UA','AU'],'28':['GU','UG']}

#===============================================================================
# find base pairs for each secondary structure
def pairing(bdstr):
    stack = []
    pairs = {}

    for i,nt in enumerate(bdstr):
        if nt == '(':
            stack.append(i)
        elif nt ==')':
            pairs[stack[-1]] = i
            pairs[i] = stack[-1]
            stack[-1:] = []
    return pairs

#===============================================================================

def getMotif(Mtype,idx,seq,bd):
    pairs = pairing(bd)
    if Mtype=='basePair':
        try:
            bp_cur = (idx,pairs[idx])
            bps = [seq[idx]+seq[pairs[idx]]]
            tps = [x for a in bps for x,y in tp.items() if a in y]
            mode = bps + tps
            index = [bp_cur[0],bp_cur[1]]
            return index,mode
        except KeyError:
            return False,False



    if Mtype=='triplet':
        try:
            bp_cur = (idx,pairs[idx])
            bp_pre = (idx-1,pairs[idx-1])
            bp_nxt = (idx+1,pairs[idx+1])
            if bp_cur[1] == bp_pre[1]-1 and bp_cur[1] == bp_nxt[1]+1:
                bps = [seq[idx]+seq[pairs[idx]],seq[idx-1]+seq[pairs[idx-1]],\
                seq[idx+1]+seq[pairs[idx+1]]] 
                tps = [x for a in bps for x,y in tp.items() if a in y]
                mode = bps + tps
                index = [bp_cur[0],bp_cur[1],bp_pre[0],bp_pre[1],bp_nxt[0],\
                bp_nxt[1]]       
                return index,mode
            else:
                return False,False
        except KeyError:
            return False,False



    if Mtype=='penta':
        try:
            bp_cur = (idx,pairs[idx])
            bp_pre = (idx-1,pairs[idx-1])
            bp_nxt = (idx+1,pairs[idx+1])
            bp_ppre = (idx-2,pairs[idx-2])
            bp_nnxt = (idx+2,pairs[idx+2])
            if bp_cur[1] == bp_pre[1]-1 and bp_cur[1] == bp_nxt[1]+1 and\
            bp_pre[1] == bp_ppre[1]-1 and bp_nxt[1] == bp_nnxt[1]+1:
                bps = [seq[idx]+seq[pairs[idx]],seq[idx-1]+seq[pairs[idx-1]],\
                seq[idx+1]+seq[pairs[idx+1]],seq[idx-2]+seq[pairs[idx-2]],\
                seq[idx+2]+seq[pairs[idx+2]]] 
                tps = [x for a in bps for x,y in tp.items() if a in y]
                mode = bps + tps
                index = [bp_cur[0],bp_cur[1],bp_pre[0],bp_pre[1],bp_nxt[0],\
                bp_nxt[1],bp_ppre[0],bp_ppre[1],bp_nnxt[0],bp_nnxt[1]]       
                return index,mode
            else:
                return False,False
        except KeyError:
            return False,False
