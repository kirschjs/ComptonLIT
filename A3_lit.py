import os, sys

import numpy as np
from scipy.optimize import fmin

from bridge import *
from rrgm_functions import *
from two_particle_functions import *
from three_particle_functions import *
from triton_width_gen import *

os.chdir(v18uixpath)

RHSofBV = {}
RHSofmJ = {}

if 'construe_fresh_helion' in cal:

    os.system('cp INLUCN_avui INLUCN')
    os.system(BINBDGpath + 'LUDW_CN.exe')
    os.system('cp INLU_avui INLU')
    os.system(BINBDGpath + 'DRLUD.exe')
    os.system('cp INOB_avui INOB')
    os.system(BINBDGpath + 'KOBER.exe')
    os.system(BINBDGpath + 'DROBER.exe')
    os.system('cp INQUA_N_V18 INQUA_N')
    os.system(BINBDGpath + 'QUAFL_N.exe')
    os.system('cp INQUA_N_UIX INQUA_N')
    os.system(BINBDGpath + 'DRQUA_AK_N.exe')
    os.system('cp INEN_avui INEN')
    os.system(BINBDGpath + 'DR2END_AK.exe')

    EBDG = get_h_ev()[0]
    np.savetxt('E0.dat', np.array([EBDG]), fmt='%12.4f')

    os.system('cp OUTPUT end_out_b && cp INEN inen_b')
    rrgm_functions.parse_ev_coeffs(infil='end_out_b')
    exit()

for file in os.listdir(litpath3He):
    if fnmatch.fnmatch(file, 'endlit*'):
        if 'dbg' in cal:
            print('removing old <*en*lit*> files.')
        os.system('rm ' + litpath3He + '/*en*lit*')
        break

for streukanal in streukas:

    os.chdir(v18uixpath)

    Jstreu = float(streukanal.split('^')[0])
    mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0, Jstreu)

    #wLIT = np.loadtxt('wLIT.dat')
    #clean_wrels = np.loadtxt('w0.dat')
    #nzfLIT = int(np.ceil(len(wLIT) / 20.0))
    #nzf0p = int(np.ceil(len(clean_wrels) / 20.0))
    #nzfTOT = nzf0p + nzfLIT

    BUECO = [cof.strip() for cof in open('COEFF')]
    EBDG = get_h_ev(ifi='end_out_b')[0]

    print(
        '(iv)    LS-scheme: B(2,%s) = %4.4f MeV [' % (boundstatekanal, EBDG),
        get_h_ev(n=4),
        ']')
    print('        dim(B_0)   = %d' % len(BUECO))
    #print('        dim(B_LIT) = %d -> %d' % (len(winiLIT), len(wLIT)))
    #print('        dim(B)     = %d' % (len(wLIT) + len(BUECO)))

    os.chdir(litpath3He)

    # assume: EINE Zerlegung pro L und S cfg. => anz int.ws < x

    lfrags = []
    sfrags = []
    for lcfg in range(len(rechtekanaele[boundstatekanal])):
        sfrags = sfrags + rechtekanaele[boundstatekanal][lcfg][1]
        for scfg in rechtekanaele[boundstatekanal][lcfg][1]:
            lfrags = lfrags + [rechtekanaele[boundstatekanal][lcfg][0]]
    for lcfg in range(len(streukanaele3He[streukanal])):
        sfrags = sfrags + streukanaele3He[streukanal][lcfg][1]
        for scfg in streukanaele3He[streukanal][lcfg][1]:
            lfrags = lfrags + [streukanaele3He[streukanal][lcfg][0]]

    lit_3inlu(mul=multipolarity, frag=lfrags)
    os.system(BINLITpath + 'luise.exe > dump')
    lit_3inob(fr=sfrags)
    os.system(BINLITpath + 'obem.exe > dump')

    intw = threebodywidths(sfrags, nPERz=6)
    relw = threebodywidths(sfrags, nPERz=4)

    lit_3inqua(anzo=11, bnd=v18uixpath + '/INQUA_N')
    lit_3inqua(intwi=intw, relwi=relw, withhead=False)
    os.system(BINLITpath + 'qual.exe')

    exit()

    leftpar = 1 if streukanal[-1] == '-' else 2

    for mM in mLmJl:
        for zerl in range(len(sfrags)):
            for bv in range(len(intw[zerl])):
                for streukanalweite in range(1, len(relwi[zerl]) + 1):
                    lit_3inen(
                        MREG='  1  0  0  0  0  0  0  0  0  1  1',
                        #                   (shifted) QBV                     nr.rw
                        KSTREU=[
                            int(streukanaeleD[streukanal][subchannel][-1] +
                                8 * nzf0p), streukanalweite
                        ],
                        JWSL=Jstreu,
                        JWSLM=mM[1],
                        MULM2=mM[0],
                        NPARL=leftpar,
                        KBND=[[2, BSRWIDX[0]], [6, BSRWIDX[1]]],
                        JWSR=J0,
                        NPARR=2,
                        EB=EBDG,
                        BUECO=BUECO,
                        NZE=anz_phot_e,
                        EK0=phot_e_0,
                        EKDIFF=phot_e_d)

                    os.system(BINLITpath + 'enemb.exe')
                    os.system('cp OUTPUT endlit%d_J%d_mJ%d-mL%d.log' %
                              (int(streukanalweite + subchannel * len(wLIT)),
                               Jstreu, mM[1], mM[0]))
                    os.system('cp INEN inenlit%d_J%d_mJ%d-mL%d.log' %
                              (int(streukanalweite + subchannel * len(wLIT)),
                               Jstreu, mM[1], mM[0]))

# basis definition
# orbital l1,l2,l12
frL = [[0, 0, 0], [1, 1, 0]]
#spin-/iso-spin
#      s12 t12  l12    T    S
# no1    1   0 even  1/2  1/2
# no2    1   0 even  1/2  3/2
# no3    0   0  odd  1/2  1/2
# no5    1   1  odd  1/2  3/2
# no6    0   1 even  1/2  1/2
frO = ['he_no1', 'he_no3']

# 2-body interaction
n3_inlu(anzO=8, fn='INLUCN', fr=frL)
os.system(BINBDGpath + 'LUDW_CN.exe')
n3_inob(frO, anzO=8, fn='INOB')
os.system(BINBDGpath + 'KOBER.exe')
repl_line('INQUA_N', 1,
          '/home/kirscher/kette_repo/sim_par/potentials/NN_pheno/AV18\n')
os.system(BINBDGpath + 'QUAFL_N.exe')

# 3-body potential
n3_inlu(anzO=8, fn='INLU', fr=frL)
os.system(BINBDGpath + 'DRLUD.exe')
n3_inob(frO, anzO=15, fn='INOB')
os.system(BINBDGpath + 'DROBER.exe')
repl_line(
    'INQUA_N', 1,
    '/home/kirscher/kette_repo/sim_par/potentials/NNN_pheno/urbana9_AK_neu\n')
os.system(BINBDGpath + 'DRQUA_AK_N.exe')

# Diagonalization
bvanz = 2
jay = 1. / 2.
rw = [1, 0]
co = [
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 0, 0, 0, 0
]
#
n3_inen_bdg(bvanz, jay, co, rw, fn='INEN', pari=0, nzop=31, tni=11)
os.system(BINBDGpath + 'DR2END_AK.exe')