import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import shutil
import re

from bridge import *


def plot_LIT_source():

    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.set_title(r'')
    ax1.set_xlabel('photon momentum [MeV]')

    [ax1.plot(photon_energy, rhs[n]) for n in range(anzcomp)]

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
