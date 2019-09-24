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


def read_uncoupled_source():

    sourceRHS = {}
    for streukanalweite in range(1, basdim + 1):

        for mM in mLmJl:
            # mM[0] = m(L) ; mM[1] = m(Jlit)

            instream = [
                line for line in open('endlitout_%d_%d-%d' % (streukanalweite,
                                                              mM[1], mM[0]))
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

                    sourceRHS[('%d' % (streukanalweite), '%d' % JLIT2,
                               '%d' % mJLIT2, '%d' % MUL2,
                               '%d' % mMUL2)] = opME
                    #print('I read r,2*(Jlit,mJlit,L,mL):', streukanalweite, JLIT2, mJLIT2, MUL2, mMUL2)

    return sourceRHS, photon_energy


# return S[ r,Jlit,m(Jlit),L, [k1,k2,...,kn] ]
def couple_source(sourceRHS):

    coupledSOURCE = {}
    for streukanalweite in range(1, basdim + 1):

        # mM[0] = m(L) ; mM[1] = m(Jlit)
        for mM in mLmJl:

            tmp = sourceRHS[('%d' % (streukanalweite),
                             '%d' % (2 * int(streukanal[0])),
                             '%d' % (2 * mM[1]), '%d' % int(2 * multipolarity),
                             '%d' % (2 * mM[0]))]
            cgtmp = CG(multipolarity, mM[0], 1, mM[1] - mM[0],
                       int(streukanal[0]), mM[1]).doit()
            #print(multipolarity, mM[0], 1, mM[1] - mM[0], int(streukanal[0]),mM[1], cgtmp)

            try:
                coupledSOURCE[('%d' % (streukanalweite),
                               '%d' % (2 * int(streukanal[0])),
                               '%d' % (2 * mM[1]),
                               '%d' % (2 * multipolarity))] += tmp * cgtmp
            except:
                coupledSOURCE[('%d' % (streukanalweite),
                               '%d' % (2 * int(streukanal[0])),
                               '%d' % (2 * mM[1]),
                               '%d' % (2 * multipolarity))] = tmp * cgtmp

    for mJ in mJlrange:

        outs = ''
        outfil = 'src_JmL_%d%d%d' % (int(streukanal[0]), mJ, multipolarity)

        for nMom in range(anz_phot_e):

            for streukanalweite in range(1, basdim + 1):
                outs += '%12.4f' % float(coupledSOURCE[(
                    '%d' % (streukanalweite), '%d' % (2 * int(streukanal[0])),
                    '%d' % (2 * mJ), '%d' % (2 * multipolarity))][nMom])

            outs += '\n'
        with open(
                av18path +
                '/LIT_SOURCE_%s%d%d' % (streukanal, mJ, multipolarity),
                'w') as outfile:
            outfile.seek(0)
            outfile.write(outs)

    return coupledSOURCE