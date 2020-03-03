import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import shutil
import re

from bridge import *
# CG(j1, m1, j2, m2, j3, m3)
from sympy.physics.quantum.cg import CG

MeVfm = 197.3161329


def read_uncoupled_source(streukanal, basisSET=''):
    # collects end OUTPUT files and stores them
    # *without* modifications

    # find mL, mJl s.t. (L,mL;Jr,mJl-mL|Jl,mJl) != 0 -------------------------
    #
    # mM[0] = m(L) ; mM[1] = m(Jlit)
    Jstreu = float(streukanal.split('^')[0])

    Jstrstr = '%s' % str(Jstreu)[:3]

    mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0, Jstreu)

    signf_overlap_bv = []
    sourceRHS = {}
    basdim = len(basisSET)

    for bv in basisSET:
        for mM in mLmJl:
            # mM[0] = m(L) ; mM[1] = m(Jlit)
            instream = [
                line
                for line in open('endlit%d-%d_J%f_mJ%f-mL%d.log' % (bv[0],
                                                                    bv[1],
                                                                    Jstreu,
                                                                    mM[1],
                                                                    mM[0]))
            ]
            #print('endlit%d-%d_J%f_mJ%f-mL%f.log' % (bv[0], bv[1], Jstreu,
            #                                         mM[1], mM[0]))
            for ln in range(len(instream)):
                if re.search('1AUSDRUCK', instream[ln]):
                    JDEUT2 = int(instream[ln + 4].split()[4])
                    JLIT2 = int(instream[ln + 4].split()[2])
                    mJLIT2 = int(instream[ln + 4].split()[5])
                    MUL2 = int(instream[ln + 4].split()[3])
                    mMUL2 = int(instream[ln + 4].split()[6])
                    photon_energy = MeVfm * np.array([
                        float(instream[ln + 4 + 3 * en].split()[1])
                        for en in range(anz_phot_e)
                    ])
                    opME = np.array([
                        float(instream[ln + 4 + 3 * en].split()[7])
                        for en in range(anz_phot_e)
                    ])

                    if (((abs(opME) > 10**-1).all()) &
                        ((abs(opME) < 10**1).all())):
                        signf_overlap_bv.append(bv)

                    sourceRHS[('%d-%d' % (bv[0],
                                          bv[1]), '%d' % JLIT2, '%d' % mJLIT2,
                               '%d' % MUL2, '%d' % mMUL2)] = opME
                    #print('I read 2*(Jlit,mJlit,L,mL):', JLIT2, mJLIT2,
                    #      MUL2, mMUL2)

    if signf_overlap_bv != []:
        with open(litpath3He + 'LITbas_red_J%s_sig.dat' % Jstrstr, 'w') as f:
            np.savetxt(f, np.unique(signf_overlap_bv, axis=0), fmt='%5d  %5d')
        f.close()

    return sourceRHS, photon_energy


# return S[ r,Jlit,m(Jlit),L, [k1,k2,...,kn] ]
def couple_source(streukanal, sourceRHS, basisSET=''):

    coupledSOURCE = {}
    basdim = len(basisSET)

    # mM[0] = m(L) ; mM[1] = m(Jlit)
    Jstreu = float(streukanal.split('^')[0])
    mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0, Jstreu)

    for bv in basisSET:
        for mJ in mJlrange:
            coupledSOURCE[('%d-%d' % (bv[0], bv[1]), '%d' % (2 * Jstreu),
                           '%d' % (2 * mJ), '%d' % (2 * multipolarity))] = 0.0

    for bv in basisSET:
        for mM in mLmJl:
            tmp = sourceRHS[('%d-%d' % (bv[0], bv[1]), '%d' % (2 * Jstreu),
                             '%d' % (2 * mM[1]), '%d' % int(2 * multipolarity),
                             '%d' % (2 * mM[0]))]
            cgtmp = CG(multipolarity, mM[0], J0, mM[1] - mM[0], Jstreu,
                       mM[1]).doit()
            coupledSOURCE[('%d-%d' % (bv[0], bv[1]), '%d' % (2 * Jstreu),
                           '%d' % (2 * mM[1]),
                           '%d' % (2 * multipolarity))] += tmp * cgtmp

    for mJ in mJlrange:

        outs = ''
        for nMom in range(anz_phot_e):

            for bv in basisSET:
                outs += '%12.4E' % float(coupledSOURCE[(
                    '%d-%d' % (bv[0], bv[1]), '%d' % (2 * Jstreu),
                    '%d' % (2 * mJ), '%d' % (2 * multipolarity))][nMom])

            outs += '\n'

        outp = litpath3He + 'LIT_SOURCE_%s_mJ' % streukanal + '{:,g}'.format(
            mJ) + '_MUL%d' % multipolarity
        if 'purge' in cal:
            if os.path.isfile(outp):
                print('removing previous <LIT_SOURCE>')
                os.system('rm ' + outp)

        with open(outp, 'w') as outfile:

            #outfile.seek(0)
            outfile.write(outs)

    return coupledSOURCE