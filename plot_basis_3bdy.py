import os
import math
import numpy as np
import matplotlib.pyplot as plt

from bridge import *


def main():

    whelionint = [float(line) for line in open(litpath3He + 'intw3he.dat')]
    whelionrel = [float(line) for line in open(litpath3He + 'relw3he.dat')]

    whelionintLIT = [
        item for line in [
            np.array(line.split(';')).astype(float).tolist()
            for line in open(litpath3He + 'intw3heLIT.dat')
        ] for item in line
    ]
    whelionrelLIT = [
        item for line in [
            np.array(line.split(';')).astype(float).tolist()
            for line in open(litpath3He + 'relw3heLIT.dat')
        ] for item in line
    ]

    rSpace = np.linspace(-20, 20, 200)

    fig = plt.figure(figsize=(10, 12), dpi=95)
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    ax1 = plt.subplot(321)
    ax1.set_title(r'Helion $\rho_1$')
    ax1.set_xlabel(r'$r$ [fm]')
    ax1.set_ylabel(r'$\phi(r)$')
    [
        ax1.plot(rSpace, [np.exp(-weite * rr**2) for rr in rSpace])
        for weite in whelionint
    ]

    ax2 = plt.subplot(322)
    ax2.set_title(r'Helion $\rho_2$')
    ax2.set_xlabel(r'$r$ [fm]')
    ax2.set_ylabel(r'$\phi(r)$')
    [
        ax2.plot(rSpace, [np.exp(-weite * rr**2) for rr in rSpace])
        for weite in whelionrel
    ]

    ax3 = plt.subplot(323)
    ax3.set_title(r'$\Psi_{LIT}\;\;\; \rho_1$')
    ax3.set_xlabel(r'$r$ [fm]')
    ax3.set_ylabel(r'$\phi(r)$')
    [
        ax3.plot(rSpace, [np.exp(-weite * rr**2) for rr in rSpace])
        for weite in whelionintLIT
    ]
    ax4 = plt.subplot(324)
    ax4.set_title(r'$\Psi_{LIT}\;\;\; \rho_1$')
    ax4.set_xlabel(r'$r$ [fm]')
    ax4.set_ylabel(r'$\phi(r)$')
    [
        ax4.plot(rSpace, [np.exp(-weite * rr**2) for rr in rSpace])
        for weite in whelionrelLIT
    ]

    ax5 = plt.subplot(313)
    ax5.set_title(r'Widths')
    ax5.set_xlabel(r'$\gamma [fm^{-2}]$')
    ax5.set_ylabel(r'')
    ax5.plot(
        whelionint,
        np.zeros(len(whelionint)),
        '.',
        label=r'$\gamma(^3He)_{\rho_1}$')
    ax5.plot(
        whelionrel,
        0.33 * np.ones(len(whelionrel)),
        '.',
        label=r'$\gamma(^3He)_{\rho_2}$')
    ax5.plot(
        whelionintLIT,
        0.66 * np.ones(len(whelionintLIT)),
        '.',
        label=r'$\gamma(\Psi_{LIT})_{\rho_1}$')
    ax5.plot(
        whelionrelLIT,
        np.ones(len(whelionrelLIT)),
        '.',
        label=r'$\gamma(\Psi_{LIT})_{\rho_2}$')
    leg = ax5.legend(loc='best')
    fig.savefig(litpath3He + 'bases.pdf')
    #plt.show()
    plt.clf()
    plt.close()


main()