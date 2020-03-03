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


def read_uncoupled_source(streukanal, basisSET=''):
    # collects end OUTPUT files and stores them
    # *without* modifications

    # find mL, mJl s.t. (L,mL;Jr,mJl-mL|Jl,mJl) != 0 -------------------------
    # ecce: Jr = Jdeuteron = 1
    # mM[0] = m(L) ; mM[1] = m(Jlit)
    mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0,
                                                  int(streukanal[0]))

    sourceRHS = {}
    basdim = len(basisSET)
    for subchannel in range(len(streukanaele[streukanal])):

        for streukanalweite in range(1, basdim + 1):

            for mM in mLmJl:
                # mM[0] = m(L) ; mM[1] = m(Jlit)

                instream = [
                    line for line in open('endlit%d_J%d_mJ%d-mL%d.log' % (
                        int(streukanalweite + subchannel * basdim),
                        int(streukanal[0]), mM[1], mM[0]))
                ]

                for ln in range(len(instream)):

                    if re.search('1AUSDRUCK', instream[ln]):
                        JDEUT2 = int(instream[ln + 3].split()[4])
                        JLIT2 = int(instream[ln + 3].split()[2])
                        mJLIT2 = int(instream[ln + 3].split()[5])
                        MUL2 = int(instream[ln + 3].split()[3])
                        mMUL2 = int(instream[ln + 3].split()[6])

                        photon_energy = np.array([
                            float(instream[ln + 3 + 2 * en].split()[1])
                            for en in range(anz_phot_e)
                        ])

                        opME = np.array([
                            float(instream[ln + 3 + 2 * en].split()[7])
                            for en in range(anz_phot_e)
                        ])

                        sourceRHS[('%d' % (
                            int(streukanalweite + subchannel * basdim)),
                                   '%d' % JLIT2, '%d' % mJLIT2, '%d' % MUL2,
                                   '%d' % mMUL2)] = opME
                        #print('I read r,2*(Jlit,mJlit,L,mL):', streukanalweite,
                        #      JLIT2, mJLIT2, MUL2, mMUL2)

    return sourceRHS, photon_energy


# return S[ r,Jlit,m(Jlit),L, [k1,k2,...,kn] ]
def couple_source(streukanal, sourceRHS, basisSET=''):

    coupledSOURCE = {}
    basdim = len(basisSET)

    # mM[0] = m(L) ; mM[1] = m(Jlit)
    mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0,
                                                  int(streukanal[0]))
    for subchannel in range(len(streukanaele[streukanal])):
        for streukanalweite in range(1, basdim + 1):
            for mJ in mJlrange:
                coupledSOURCE[('%d' %
                               (int(streukanalweite + subchannel * basdim)),
                               '%d' % (2 * int(streukanal[0])), '%d' %
                               (2 * mJ), '%d' % (2 * multipolarity))] = 0.0

    for subchannel in range(len(streukanaele[streukanal])):

        for streukanalweite in range(1, basdim + 1):

            for mM in mLmJl:

                tmp = sourceRHS[('%d' %
                                 (int(streukanalweite + subchannel * basdim)),
                                 '%d' % (2 * int(streukanal[0])), '%d' %
                                 (2 * mM[1]), '%d' % int(2 * multipolarity),
                                 '%d' % (2 * mM[0]))]
                cgtmp = CG(multipolarity, mM[0], 1, mM[1] - mM[0],
                           int(streukanal[0]), mM[1]).doit()

                #print(multipolarity, mM[0], 1, mM[1] - mM[0], int(streukanal[0]), mM[1], cgtmp)

                coupledSOURCE[('%d' %
                               (int(streukanalweite + subchannel * basdim)),
                               '%d' % (2 * int(streukanal[0])),
                               '%d' % (2 * mM[1]),
                               '%d' % (2 * multipolarity))] += tmp * cgtmp

    for mJ in mJlrange:

        outs = ''
        for nMom in range(anz_phot_e):

            for subchannel in range(len(streukanaele[streukanal])):

                for streukanalweite in range(1, basdim + 1):

                    outs += '%12.4E' % float(coupledSOURCE[(
                        '%d' % (int(streukanalweite + subchannel * basdim)),
                        '%d' % (2 * int(streukanal[0])), '%d' % (2 * mJ),
                        '%d' % (2 * multipolarity))][nMom])

            outs += '\n'

        outp = av18path + '/LIT_SOURCE_%s%d%d' % (streukanal, mJ,
                                                  multipolarity)
        if 'purge' in cal:
            if os.path.isfile(outp):
                print('removing previous <LIT_SOURCE>')
                os.system('rm ' + outp)

        with open(outp, 'w') as outfile:

            #outfile.seek(0)
            outfile.write(outs)

    return coupledSOURCE