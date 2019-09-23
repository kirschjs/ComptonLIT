import os
import numpy as np
import sympy as sy
# CG(j1, m1, j2, m2, j3, m3)
from sympy.physics.quantum.cg import CG

from parameters_and_constants import *
from rrgm_functions import *
from two_particle_functions import *

pathbase = '/home_th/kirscher/kette_repo/source/ComptonLIT'

av18path = pathbase + '/av18_deuteron'
litpath = pathbase + '/test/'
BINBDGpath = pathbase + '/src_nucl/'
BINLITpath = pathbase + '/src_elma/'

mpii = '137'

cal = ['bdg', 'QUA']
cal = ['bdg', 'QUA', 'reduce']
dbg = 2

pot_scale = 1.

anzs = 1
v_i = 1.0
v_e = 1.0
over_space = [15., 0.1, 0.0001]
eps_space = np.linspace(v_i, v_e, anzs)

anze = 17

multipolarity = 1
mLrange = np.arange(-multipolarity, multipolarity + 1)

anz_phot_e = 100
phot_e_0 = 0.1  #  enems_e converts to fm^-1, but HERE the value is in MeV
phot_e_d = 10.  #

addw = 9
scale = 1.1
wini = w120
sizeFrag = len(wini) - 1
basdim = addw + len(wini)

boundstatekanal = 'np-3SD1'
streukanal = '0-'
mJlrange = np.arange(-int(streukanal[0]), int(streukanal[0]) + 1)

# find mL, mJl s.t. (L,mL;Jr,mJl-mL|Jl,mJl) != 0
# ecce: Jr = Jdeuteron = 1
mLmJl = []
for mM in np.array(np.meshgrid(mLrange, mJlrange)).T.reshape(-1, 2):
    clg = CG(multipolarity, mM[0], 1, mM[1] - mM[0], int(streukanal[0]), mM[1])
    if (clg.doit() == 0):
        continue
    mLmJl.append(mM)
print(mLmJl)
exit()
streukanalweiten = range(1, len(wini) + addw + 1)
#wdim = 20
#loc = 1e-4
#wini = abs(np.random.laplace(loc, scale, wdim))

maxCoef = 10000
minCoef = 1200
ncycl = 30
maxDiff = 0.001
delPcyc = 1