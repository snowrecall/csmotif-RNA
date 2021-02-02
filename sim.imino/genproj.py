#!/usr/bin/python

from nmrglue import *
from mplot import *

#xl = [14.5, 10.3]
#yl = [0, 12000000]

pagesize('l')

dic,data = sparky.read('./sim.ucsf')
uc = sparky.make_uc(dic, data, dim=1)
plot(uc.ppm_scale(), data.max(0),'k-')
pyl.gca().invert_xaxis()
#xlim(xl)
#ylim(yl)

saveplot('proj.pdf')
show()
