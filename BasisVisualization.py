from bridge import *
from three_particle_functions import *
from triton_width_gen import *
import operator

import matplotlib.pyplot as plt


def visbas(basispath, widthpath, exepath, Jstrstr='0.5'):

    curd = os.getcwd()
    os.chdir(basispath)

    litbas_full = np.loadtxt(
        basispath + 'LITbas_full_J%s.dat' % Jstrstr).astype(int)
    try:
        litbas_red = np.loadtxt(
            basispath + 'LITbas_red_J%s.dat' % Jstrstr).astype(int)
    except:
        litbas_red = litbas_full

    basvs = {}
    for bv in litbas_red:
        try:
            basvs[bv[0]].append(bv[1])
        except:
            basvs[bv[0]] = [bv[1]]
    sbas = sorted(basvs.items(), key=operator.itemgetter(0))

    n3_inen_bdg(
        sbas, float(Jstrstr), costr, fn='INEN', pari=0, nzop=31, tni=11)

    intwLIT = [
        np.array(ln.split(';')).astype(float).tolist()
        for ln in open(widthpath + 'intw3heLIT_J%s.dat' % Jstrstr)
    ]

    relwLIT = [
        np.array(ln.split(';')).astype(float).tolist()
        for ln in open(widthpath + 'relw3heLIT_J%s.dat' % Jstrstr)
    ]

    iws_full = []
    rws_full = []

    iws_red = []
    rws_red = []

    numered_widths = []
    #print(' BV REL          wi          wr')

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

                numered_widths.append([bv[0], bv[1], iws_red[-1], rws_red[-1]])
                break

    numered_widths = np.array(numered_widths)

    numered_widths = numered_widths[np.lexsort(([
        numered_widths[:, i] for i in range(numered_widths.shape[1] - 2,
                                            numered_widths.shape[1] - 1, +1)
    ]))]

    with open(v18uixpath + 'LITbas_red_J%s_doc.dat' % Jstrstr, 'w') as f:
        np.savetxt(f, numered_widths, fmt='%5d  %5d  %12.8f  %12.8f')
    f.close()

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

    plt.title(
        r'$\gamma_1<\gamma_2\Rightarrow$ 1 broader than 2   (J=%s)' % Jstrstr)

    fig.savefig(basispath + 'WidthXY_J%s.pdf' % Jstrstr)

    os.chdir(curd)


#Jstreu = float(streukas[-1].split('^')[0])
#Jstreustring = '%s' % str(Jstreu)[:3]
#
#os.chdir(v18uixpath)
#
#os.system('cp QUAOUT_J%s QUAOUT' % Jstreustring)
#os.system('cp DRQUAOUT_J%s DRQUAOUT' % Jstreustring)
#
#visbas(
#    basispath=v18uixpath,
#    widthpath=litpath3He,
#    exepath=BINBDGpath,
#    Jstrstr=Jstreustring)
#
#os.system(BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')