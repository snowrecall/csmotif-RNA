#! /usr/bin/python

import nmrglue as ng
import numpy as np
from os import path
from common.base import divide
from sparky import *

NU = [('H1','N1'), ('H3','N3')]

# read in the CS data
lines = open('../imino.tab').readlines()
lines = filter(lambda x: x[:1]!='#', lines)
# build lookup table
data = {}
blks = divide(lines, lambda x: x.split()[1])
for blk in blks:
    fds = blk[0].split()
    resn,resi = fds[0],int(fds[1])
    CS = {}
    for line in blk:
        fds = line.split()
        resn,resi,nu,cs = fds[0],int(fds[1]),fds[2],float(fds[3])
        CS[nu] = cs
    data[resi] = [resn, CS]
# find peaks
pks = []
resis = data.keys()
resis.sort()
for resi in resis:
    resn,CS = data[resi]
    for nu1,nu2 in NU:
        if nu1 in CS and nu2 in CS:
            pks.append([resn+`resi`,CS[nu1],CS[nu2]])
print 'simulated %d imino peaks: '%len(pks) + ' '.join(map(lambda x: x[0], pks))

npeaks = len(pks)
csH = np.array([item[1] for item in pks])
csN = np.array([item[2] for item in pks])

#sw = csH.max() - csH.min()
#carH = (csH.min()+sw/2.0) * 600.0
#swH = sw * 1.2 * 600.0
#sw = csN.max() - csN.min()
#carN = (csN.min()+sw/2.0) * 60.8
#swN = sw * 1.2 * 60.8
maxcsH,mincsH = 14.5, 10.
maxcsN,mincsN = 165., 140.
sw = maxcsH - mincsH
carH = (mincsH+sw/2.0) * 600.0
swH = sw * 1.2 * 600.0
sw = maxcsN - mincsN
carN = (mincsN+sw/2.0) * 60.8
swN = sw * 1.2 * 60.8

# create a sparky dictionary
# A dictionary from a existing Sparky ucsf file can be found using:
# ng.sparky.guess_udic(*ng.sparky.read('filename.ucsf'))
udic = {
    'ndim': 2,
    0: {'car': carN,
        'complex': False,
        'encoding': 'states',
        'freq': True,
        'label': '15N',
        'obs': 60.8,
        'size': 512,
        'sw': swN,
        'time': False},
    1: {'car': carH,
        'complex': False,
        'encoding': 'direct',
        'freq': True,
        'label': '1H',
        'obs': 600.0,
        'size': 1024,
        'sw': swH,
        'time': False}
}

dic = ng.sparky.create_dic(udic)
data = np.empty((512, 1024), dtype='float32')

# convert the peak list from PPM to points
uc_15N = ng.sparky.make_uc(dic, None, 0)
uc_1H = ng.sparky.make_uc(dic, None, 1)

lw_15N = 6.0    # 15N dimension linewidth in points
lw_1H = 12.0     # 1H dimension linewidth in points

params = []
for ass,ppm_1H, ppm_15N in pks:
    pts_15N = uc_15N.f(ppm_15N, 'ppm')
    pts_1H = uc_1H.f(ppm_1H, 'ppm')
    params.append([(pts_15N, lw_15N), (pts_1H, lw_1H)])

# simulate the spectrum
shape = (512, 1024)      # size should match the dictionary size
lineshapes = ('g', 'g')  # gaussian in both dimensions
amps = [100.0] * npeaks
data = ng.linesh.sim_NDregion(shape, lineshapes, params, amps)

# save the spectrum
ng.sparky.write("sim.ucsf", dic, data.astype('float32'), overwrite=True)

# generate a save file if possible
if path.isfile('tpl.save'):
    sv = SVFile('tpl.save')
    i = 1
    for ass,csx,csy in pks:
        pk = Node()
        pk.id = i
        pk.label = ass
        pk.rs = [[ass[0], ass[1:]], ['?', '?']]
        pk.labelpos = [csy, csx]
        pk.labelxy = [[csy, csx-0.035], [csx+0.16, csy-1.0]]
        pk.pos = [csy, csx]
        pk.height = [0.0, 100]
        sv.pks.append(pk)
        i += 1
    sv.write('sim.save')
