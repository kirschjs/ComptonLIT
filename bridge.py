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
cal = ['lit']
cal = ['lit-plot']
cal = ['bdg', 'QUA', 'reduce', 'lit', 'lit-plot']

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

# deuteron/initial-state basis -------------------------------------------
basisdim0 = 35

laplace_loc, laplace_scale = 1., .4
wLAPLACE = np.sort(
    np.abs(np.random.laplace(laplace_loc, laplace_scale, basisdim0)))
wini0 = wLAPLACE

addw = 8
addwt = 'middle'
scale = 1.
min_spacing = 0.6

print('initial width set : ', wini0)
rw0 = wid_gen(
    add=addw, addtype=addwt, w0=wini0, ths=[1e-5, 2e2, 0.2], sca=scale)
rw0 = sparsify(rw0, min_spacing)
print('extended width set: ', rw0)
nzf0 = int(np.ceil(len(rw0) / 20.0))
print('NZF = ', nzf0)

#LIT basis ---------------------------------------------------------------

basdim = addw + len(wini)

basisdimLIT = 35
winiLIT = np.geomspace(
    start=0.1, stop=100, num=basisdimLIT, endpoint=True, dtype=None)

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

streukanalweiten = range(1, len(wini) + addw + 1)
#wdim = 20
#loc = 1e-4
#wini = abs(np.random.laplace(loc, scale, wdim))

maxCoef = 10000
minCoef = 1200
ncycl = 30
maxDiff = 0.001
delPcyc = 1