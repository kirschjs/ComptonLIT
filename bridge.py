import numpy as np

from parameters_and_constants import *
from rrgm_functions import *
from two_particle_functions import *

from pathlib import Path
homedir = str(Path.home())

av18path = homedir + '/kette_repo/sim_par/nucleus/2n/LIT_prep/beyond_siegert/av18_deuteron'
BINBDGpath = homedir + '/kette_repo/source/seriell/untainted/'
BINLITpath = homedir + '/kette_repo/source/seriell/lit_rrgm_elmag/new_noentw/'

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
anz_phot_e = 100
phot_e_0 = 0.1  #  enems_e converts to fm^-1, but HERE the value is in MeV
phot_e_d = 10.  #

addw = 9
scale = 1.1
wini = w120
sizeFrag = len(wini) - 1

boundstatekanal = 'np-3SD1'
streukanal = '0-'
streukanalweiten = range(1, len(wini) + addw + 1)
#wdim = 20
#loc = 1e-4
#wini = abs(np.random.laplace(loc, scale, wdim))

maxCoef = 10000
minCoef = 1200
ncycl = 30
maxDiff = 0.001
delPcyc = 1