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
pots = 'AV18'

# convention: bound-state-expanding BVs: (1-8), i.e., 8 states per rw set => nzf0*8
rechtekanaele = {
    # DEUTERON
    'np-3SD1': [[1, 1, 2, 2], [1, 1, 2, 6]],
}

# ECCE! bvnr for LIT basis states must be shifted after purge! (A2_lit)
streukanaele = {
    #       s1 s2 S  bvnr
    '0+': [[1, 1, 0, 1]],  #     1S0
    '0-': [[1, 1, 2, 4]],  #     3P0
    '1+': [[1, 1, 2, 2], [1, 1, 2, 6]],  # 3S1 3D1
    #'1-': [[1, 1, 0, 3], [1, 1, 2, 4]],  # 1P1 3P1
    '1-': [[1, 1, 2, 4]],  # 3P1
    '2+': [[1, 1, 0, 5], [1, 1, 2, 6]],  # 1D2 3D2
    #'2-': [[1, 1, 2, 4], [1, 1, 2, 8]],  # 3P2 3F2
    '2-': [[1, 1, 2, 4]],  # 3P2
    '3+': [[1, 1, 2, 6]],  #     3D3
    '3-': [[1, 1, 0, 7], [1, 1, 2, 8]],  # 1F3 3F3
    '4-': [[1, 1, 2, 8]],  #     3F4
}

cal = ['bdg', 'QUA', 'lit', 'lit-plot']

multipolarity = 1
mLrange = np.arange(-multipolarity, multipolarity + 1)

anz_phot_e = 100
phot_e_0 = 0.2  #  enems_e converts to fm^-1, but HERE the value is in MeV
phot_e_d = 10.  #

# deuteron/initial-state basis -------------------------------------------
basisdim0 = 35

laplace_loc, laplace_scale = 1., .4
wLAPLACE = np.sort(
    np.abs(np.random.laplace(laplace_loc, laplace_scale, basisdim0)))
wini0 = wLAPLACE[::-1]

addw = 8
addwt = 'middle'
scale = 1.
min_spacing = 0.2

rw0 = wid_gen(
    add=addw, addtype=addwt, w0=wini0, ths=[1e-5, 2e2, 0.2], sca=scale)
rw0 = sparsify(rw0, min_spacing)

nzf0 = int(np.ceil(len(rw0) / 20.0))

#LIT basis ---------------------------------------------------------------

basisdimLIT = 10
#winiLITa = np.geomspace(
#    start=1.0001, stop=20, num=basisdimLIT, endpoint=True, dtype=None)
laplace_loc, laplace_scale = 1., .9
winiLIT = np.sort(
    np.abs(np.random.laplace(laplace_loc, laplace_scale, basisdimLIT)))

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

streukanalweiten = range(1, len(winiLIT) + 1)
#wdim = 20
#loc = 1e-4
#wini = abs(np.random.laplace(loc, scale, wdim))

maxCoef = 10000
minCoef = 1200
ncycl = 30
maxDiff = 0.001
delPcyc = 1