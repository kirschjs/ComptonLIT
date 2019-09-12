import os, sys

import numpy as np
from scipy.optimize import fmin

from bridge import *
from rrgm_functions import *
from two_particle_functions import *

os.chdir(av18path)

pots = 'AV18'

rechtekanaele = {
    # DEUTERON
    'np-3SD1': [[1, 1, 2, 2], [1, 1, 2, 6]],
}
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

os.system(BINBDGpath + 'LUDW_CN.exe')
os.system(BINBDGpath + 'KOBER.exe')

phasSeum = [np.zeros(anze) for n in range(anzs)]

tmp1 = []
tmp2 = []
for n in range(len(eps_space)):

    if ('QUA' in cal):
        rw = wid_gen(add=addw, w0=wini, ths=[1e-5, 2e2, 0.2], sca=scale)
        h2_inqua(rw, pots)
        os.system(BINBDGpath + 'QUAFL_N.exe')
        costr = ''
        for nn in range(14):
            costr += '%12.7f' % 1.0 if (nn != 6) else '%12.7f\n' % 1.0
        h2_inen_bs(
            relw=rw,
            costr=costr,
            j=int(boundstatekanal[-1]),
            ch=rechtekanaele[boundstatekanal])
    else:
        rw = get_quaf_width_set()

    if (('bdg' in cal) | ('reduce' in cal)):
        os.system(BINBDGpath + 'DR2END_AK.exe')
        EBDG = get_h_ev()[0]
        if dbg:
            print(
                'LS-scheme: B(2,%s,eps=%2.2f) = %4.4f MeV [' %
                (boundstatekanal, eps_space[n], EBDG),
                get_h_ev(n=4),
                ']')
        rrgm_functions.parse_ev_coeffs()
        os.system('cp OUTPUT end_out_b && cp INEN inen_b')

    if 'reduce' in cal:
        reduce_2n(
            w2rels=rw,
            ch=boundstatekanal,
            size2=sizeFrag,
            ncycl=ncycl,
            maxd=maxDiff,
            minc2=minCoef,
            maxc2=maxCoef,
            execut=homedir +
            '/kette_repo/source/seriell/untainted/DR2END_AK.exe')
        if dbg:
            print('-- reduced B(2,%s) = %4.4f MeV' % (boundstatekanal,
                                                      get_h_ev()[0]))
        os.system('cp OUTPUT end_out_b && cp INEN inen_b')

rrgm_functions.parse_ev_coeffs()
BUECO = [cof.strip() for cof in open('COEFF')]
BSRWIDX = get_bsv_rw_idx(chs=2)

litpath = homedir + '/kette_repo/source/seriell/lit_rrgm_elmag/new_noentw/test/'
os.chdir(litpath)

lit_inlu(mul=multipolarity)
os.system(BINLITpath + 'luise.exe')
lit_inob()
os.system(BINLITpath + 'obem.exe')
lit_inqua(relw=rw)
os.system(BINLITpath + 'qual.exe')

leftpar = 1 if streukanal[1] == '-' else 2

rwn = 0
for streukanalweite in range(1, len(rw) - 1):
    rwn += 1
    lit_inen(
        MREG='  1  0  0  0  0  0  0  0  0  1  1',
        #                                  QBV  nr.rw
        KSTREU=[streukanaele[streukanal][0][-1], streukanalweite],
        JWSL=int(streukanal[0]),
        NPARL=leftpar,
        KBND=[[2, BSRWIDX[0]], [5, BSRWIDX[1]]],
        JWSR=1,
        NPARR=2,
        rw=rw,
        EB=EBDG,
        BUECO=BUECO,
        NZE=anz_phot_e,
        EK0=phot_e_0,
        EKDIFF=phot_e_d)

    os.system(BINLITpath + 'enemb.exe')
    os.system('cp OUTPUT endlitout_%d' % rwn)

os.system(
    'python /home_th/kirscher/kette_repo/source/seriell/lit_rrgm_elmag/new_noentw/lit_plots.py'
)

print('calculating norm/ham in scattering-channel basis')
os.chdir(av18path)
h2_inen_bs(
    relw=rw, costr=costr, j=int(streukanal[0]), ch=streukanaele[streukanal])
os.system(BINBDGpath + 'DR2END_AK.exe')
os.system(
    'cp %s/MATOUT %s/norm-ham-litME-%s' % (av18path, av18path, streukanal))