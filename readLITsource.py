import os
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import shutil
import re

from bridge import *


def read_uncoupled_source():

    for streukanalweite in range(1, basdim + 1):

        for mM in mLmJl:

            instream = [
                line for line in open('endlitout_%d_%d-%d' % (streukanalweite,
                                                              mM[1], mM[0]))
            ]
            for ln in range(len(instream)):
                print(instream[ln])
                if re.search('1AUSDRUCK', instream[ln]):
                    JDEUT = int(instream[ln + 3].split()[4])
                    JLIT = int(instream[ln + 3].split()[2])
                    mJLIT = int(instream[ln + 3].split()[5])
                    MUL = int(instream[ln + 3].split()[3])
                    mMUL = int(instream[ln + 3].split()[6])

                    photon_energy = [
                        float(instream[ln + 3 + 2 * en].split()[1])
                        for en in range(anz_phot_e)
                    ]
                    opME = [
                        float(instream[ln + 3 + 2 * en].split()[7])
                        for en in range(anz_phot_e)
                    ]
            return JLIT, mJLIT, MUL, mMUL, JDEUT, photon_energy, opME