from bridge import *
from three_particle_functions import *
from triton_width_gen import *
import operator

import matplotlib.pyplot as plt


def visbas(basispath, widthpath, exepath):

    os.chdir(basispath)

    litbas_full = np.loadtxt(basispath + 'LITbas_full.dat').astype(int)
    litbas_red = np.loadtxt(basispath + 'LITbas_red.dat').astype(int)

    basvs = {}
    for bv in litbas_red:
        try:
            basvs[bv[0]].append(bv[1])
        except:
            basvs[bv[0]] = [bv[1]]
    sbas = sorted(basvs.items(), key=operator.itemgetter(0))

    n3_inen_bdg(sbas, 1. / 2., costr, fn='INEN', pari=0, nzop=31, tni=11)

    intwLIT = [
        np.array(ln.split(';')).astype(float).tolist()
        for ln in open(widthpath + 'intw3heLIT.dat')
    ]

    relwLIT = [
        np.array(ln.split(';')).astype(float).tolist()
        for ln in open(widthpath + 'relw3heLIT.dat')
    ]

    iws_full = []
    rws_full = []

    iws_red = []
    rws_red = []

    print(' BV REL          wi          wr')

    for bv in litbas_full:
        for fr in range(len(intwLIT)):
            if bv[0] <= sum([len(ws) for ws in intwLIT[:fr + 1]]):
                iw = bv[0] - sum([len(ws) for ws in intwLIT[:fr]])
                iws_full.append(float(intwLIT[fr][iw - 1]))
                rws_full.append(float(relwLIT[fr][bv[1] - 1]))
                break

    for bv in litbas_red:
        for fr in range(len(intwLIT)):
            if bv[0] <= sum([len(ws) for ws in intwLIT[:fr + 1]]):
                iw = bv[0] - sum([len(ws) for ws in intwLIT[:fr]])
                iws_red.append(float(intwLIT[fr][iw - 1]))
                rws_red.append(float(relwLIT[fr][bv[1] - 1]))
                print('%3d%3d%12.4f%12.4f' % (bv[0], bv[1], iws_red[-1],
                                              rws_red[-1]))
                break

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.plot(
        iws_full,
        rws_full,
        'o',
        alpha=0.3,
        color='blue',
        label=r'full basis set',
        marker='.',
        markersize=2)
    ax.plot(
        iws_red,
        rws_red,
        'o',
        alpha=0.8,
        marker='.',
        markeredgecolor='red',
        markersize=4,
        color='k',
        label=r'reduced basis set')

    plt.xlabel(r'$\gamma(\rho_1)\;\;\; [fm^{-2}]$')
    plt.ylabel(r'$\gamma(\rho_2)\;\;\; [fm^{-2}]$')

    plt.legend(loc='best')

    plt.title(r'$\gamma_1<\gamma_2\Rightarrow$ 1 broader than 2')

    fig.savefig(basispath + 'WidthXY.pdf')
    os.system(exepath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')


visbas(basispath=v18uixpath, widthpath=litpath3He, exepath=BINBDGpath)