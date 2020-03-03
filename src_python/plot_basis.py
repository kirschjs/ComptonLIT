import os
import math
import numpy as np
import matplotlib.pyplot as plt


def main():
    home = os.getenv("HOME") + '/kette_repo/ComptonLIT/av18_deuteron/'

    wDeuteron = [float(line) for line in open(home + 'w0.dat')]
    wLIT = [float(line) for line in open(home + 'wLIT.dat')]

    rSpace = np.linspace(-10, 10, 200)

    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4, wspace=0.4)

    ax1 = plt.subplot(221)
    ax1.set_title(r'Deuteron')
    ax1.set_xlabel(r'$r$ [fm]')
    ax1.set_ylabel(r'$\phi(r)$')
    [
        ax1.plot(rSpace, [np.exp(-weite * rr**2) for rr in rSpace])
        for weite in wDeuteron
    ]
    #ax.text(0.5, 0.5, str((2, 3, i)),
    #       fontsize=18, ha='center')
    ax2 = plt.subplot(222)
    ax2.set_title(r'LIT')
    ax2.set_xlabel(r'$r$ [fm]')
    ax2.set_ylabel(r'$\phi(r)$')
    [
        ax2.plot(rSpace, [np.exp(-weite * rr**2) for rr in rSpace])
        for weite in wLIT
    ]

    ax3 = plt.subplot(212)
    ax3.set_title(r'Widths')
    ax3.set_xlabel(r'$\gamma [fm^{-2}]$')
    ax3.set_ylabel(r'')
    ax3.plot(wLIT, np.zeros(len(wLIT)), '.', label=r'$\gamma_{LIT}$')
    ax3.plot(
        wDeuteron, np.ones(len(wDeuteron)), '.', label=r'$\gamma_{Deuteron}$')
    leg = ax3.legend(loc='best')
    fig.savefig(home + 'bases.pdf')
    #plt.show()
    plt.clf()
    plt.close()