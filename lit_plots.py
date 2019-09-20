import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import shutil
import re

from bridge import *


def plot_LIT_source(mM=mLmJl[0]):
    pltfile = 'OUTPUT'

    anzcomp = len(streukanalweiten)
    endfile = 'endlitout_(%d->%d)_%d-%d' % (1, anzcomp, mM[1], mM[0])
    print('(iiia)  plotting %s' % endfile)

    outputend = [[
        line
        for line in open(litpath + 'endlitout_%d_%d-%d' % (cmp, mM[1], mM[0]))
    ] for cmp in range(1, anzcomp + 1)]

    mul = int([ln for ln in open(litpath + 'INLU')][2].split()[1])
    anzmom = int([ln for ln in open(litpath + 'INEN')][6])

    if ((mul != multipolarity) | (anzmom != anz_phot_e)):
        print(mul, multipolarity, anzmom, anz_phot_e)
        print('Data inconsistency! exiting...')
        exit()

    photon_energy = np.empty(anzmom)
    rhs = np.empty([anzcomp, anzmom])
    rhsb = np.empty([anzcomp, anzmom])

    for ln in range(len(outputend[0])):
        if re.search('1AUSDRUCK', outputend[0][ln]):
            photon_energy = [
                float(outputend[0][ln + 3 + 2 * en].split()[1])
                for en in range(anzmom)
            ]

    for n in range(len(outputend)):
        for ln in range(len(outputend[n])):
            if re.search('1AUSDRUCK', outputend[n][ln]):
                for en in range(anzmom):
                    rhs[n][en] = float(
                        outputend[n][ln + 3 + 2 * en].split()[-2])
                    rhsb[n][en] = float(
                        outputend[n][ln + 3 + 2 * en + 1].split()[-2])

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.set_title(r'')
    ax1.set_xlabel('photon momentum [MeV]')

    [ax1.plot(photon_energy, rhs[n]) for n in range(anzcomp)]
    [ax2.plot(photon_energy, rhsb[n]) for n in range(anzcomp)]

    plt.show()


#outstr_head = '# k_photon [MeV]'
#for bvn in range(anzcomp):
#    outstr_head += '%13s' % str(bvn)
#
#outstr_head += '\n'
#
#outstr = ''
#
#for en in range(anzmom):
#    outstr += '%15s' % str(photon_energy[en])
#    for bvn in range(anzcomp):
#        outstr += '%12.4e ' % (rhsb[bvn][en])
#    outstr += '\n'
#
#with open(av18path + '/LIT_SOURCE-%s' % streukanal, 'w') as outfile:
#    outfile.seek(0)
#    outfile.write(outstr)
