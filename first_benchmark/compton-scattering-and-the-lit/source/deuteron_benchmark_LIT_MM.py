import os, sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate
from math import factorial as fac
from scipy.special import comb

from sympy.physics.wigner import wigner_3j

# functions: clg(5I), f6j(6I), f9j(9I), s6j(6I)
# f6j performs triangle test before invoking s6j
# hmh's functions take INTs => hand over 2*S
# test: 6j(1,0,1\\1,0,1)=1/3
#import wign

MeVfm = 197.3161329

# read polarizabilities from Winfried's files for
# dipole-approximated compton-scattering on a H2
polarizability = {}
kmin = 0
kmax = 1e24
ktest = 40.97
nbr_energies = 100

print('k = %2.1f MeV' % ktest)
print('       Re           Im')

for jay in range(0, 3):
    print('P_J%d = ' % (jay), end='')

    for comp in ['Re', 'Im']:
        #joffset = [1.25 * 2.11016415E-03, 0,
        #           -0.945 * 3.70857771E-04] if comp == 'Re' else [0., 0., 0.]
        file_str = './LIT_results/lit-pol-prc/' + comp + '_P_J' + str(jay)
        if jay == 1:
            file_str += '_mod'
        tmp = np.array([
            np.array(line.strip().split()).astype(float)
            for line in open(file_str)
        ])
        kmin = max(kmin, min(tmp[:, 0]))
        kmax = min(kmax, max(tmp[:, 0]))
        polarizability[jay, comp] = scipy.interpolate.interp1d(
            tmp[:, 0], tmp[:, 1], kind='cubic')
        print('%+-5.7f   ' % polarizability[jay, comp](ktest), end='')
    print('')


def w3j(a, b, c, d, e, f):

    ww = float(wigner_3j(1 * a, 1 * b, 1 * c, 1 * d, 1 * e, 1 * f))

    return ww
    #return (clebgor(a, b, c, d, e,(d + e)) / np.sqrt(2 * c + 1)) * (-1)**(a - b - f)


def rotMat(j, mp, m, beta):

    #if ((abs(m) > abs(j)) | (abs(mp) > abs(j))):
    #    return 0
    try:
        c1 = np.sqrt(fac(j + mp) * fac(j - mp) / (fac(j + m) * fac(j - m)))
    except:
        c1 = 0
    summand = 0

    for sigma in range(j + abs(m) + 1):

        combi = comb(j + m, j - mp - sigma) * comb(j - m, sigma)
        if combi != 0:
            summand += combi * (-1)**(j - sigma - mp) * np.cos(
                0.5 * beta)**(2 * sigma + m + mp) * np.sin(
                    0.5 * beta)**(2 * j - 2 * sigma - m - mp)

    return c1 * summand


def plot_pol(jj=0):

    xl = {
        0: [[-3.4e-3, -0.9e-3], [0, 1.8e-3]],
        1: [[-1.e-4, 5.e-4], [-9.2e-5, 2.e-5]],
        2: [[0., 4.e-4], [-2.5e-4, 0.]]
    }

    fig = plt.figure()
    ax1 = fig.add_subplot(121)

    ax1.set_xlabel(r'$k\;\;[MeV]$', fontsize=16)

    xx = np.linspace(kmin, kmax, nbr_energies)
    yy = polarizability[jj, 'Im'](xx)

    plt.ylim(xl[jj][1])
    ax1.set_ylabel(r'$Im(P_{if,J=%d})\;\;[fm]$' % (jj), fontsize=14)

    ax1.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    ax1.tick_params(axis='y', colors='b')

    ax1.plot(xx, yy, '-', lw=2, alpha=1., label=r'%s' % ('Im'), color='b')

    ax1 = fig.add_subplot(122)
    yy = polarizability[jj, 'Re'](xx)

    plt.ylim(xl[jj][0])
    ax1.set_ylabel(r'$Re(P_{if,J=%d})\;\;[fm]$' % (jj), fontsize=14)
    ax1.ticklabel_format(
        useOffset=False, axis='y', style='sci', scilimits=(-2, 1))

    ax1.plot(xx, yy, linestyle='--', lw=1.5, alpha=1., color='r')

    plt.show()


def calc_amplitude(kk=40,
                   lambd=1,
                   lambdap=-1,
                   theta=np.pi / 4,
                   Mi=0,
                   Mf=0,
                   jset=[0, 1, 2],
                   jphase=[],
                   phase=[1, 1, 1],
                   pola=[]):

    summand = 0
    for jj in jset:

        polari = phase[jj] * (-1)**(lambdap - Mi) * (
            polarizability[jj, 'Re'](kk) + 1j * polarizability[jj, 'Im']
            (kk)) if pola == [] else pola[jj]
        summand += (2 * jj + 1) * w3j(1, jj, 1, -Mf, Mf - Mi, Mi) * w3j(
            1, 1, jj, lambd, Mf - Mi - lambd, Mi - Mf) * polari * rotMat(
                1, (Mf - Mi - lambd), -lambdap, -theta)

    return summand


def plot_ampl(la,
              lap,
              mi,
              mf,
              th,
              jselection=[0, 1, 2],
              jphas=[],
              ph=[1, 1, 1]):

    fig, ax1 = plt.subplots(
        num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')

    ax1.set_xlabel(r'$k\;\;[MeV]$', fontsize=16)

    krange = np.linspace(kmin, kmax, nbr_energies)

    imagp = np.imag(
        calc_amplitude(
            kk=krange,
            lambd=la,
            lambdap=lap,
            theta=th,
            Mi=mi,
            Mf=mf,
            jset=jselection,
            jphase=jphas,
            phase=ph))

    ax1.set_ylabel(
        r'$%s(T^{m_i=%d,m_f=%d}_{%d})\;\;[fm]$' % ('Im', mi, mf, lap),
        fontsize=14)

    ax1.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    ax1.tick_params(axis='y', colors='b')
    ax1.plot(imagp[:, 0], imagp[:, 1], '-', lw=2, alpha=.5, color='b')

    realp = np.real(
        calc_amplitude(
            kk=krange,
            lambd=la,
            lambdap=lap,
            theta=th,
            Mi=mi,
            Mf=mf,
            jset=jselection,
            jphase=jphas,
            phase=ph))

    ax2 = ax1.twinx()
    ax2.set_ylabel(
        r'$%s(T^{m_i=%d,m_f=%d}_{%d})\;\;[fm]$' % ('Re', mi, mf, lap),
        fontsize=14)

    ax2.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    ax2.tick_params(axis='y', colors='r')

    ax2.plot(realp[:, 0], realp[:, 1], '-', lw=2, alpha=.5, color='r')

    plt.title(r'$\theta=%2.2f^{\circ}$' % (th * 180 / np.pi))
    plt.show()


def plot_ampl_of_theta(kk,
                       indizes=[[1, 1, 1, 1]],
                       nbr_angles=20,
                       jselection=[0, 1, 2],
                       jphas=[],
                       ph=[1, 1, 1],
                       pola=[]):

    scattering_angle_range = np.linspace(0.0, np.pi, nbr_angles)

    ampl_of_theta = []

    for ampi in indizes:
        tmp = []
        for thet in scattering_angle_range:
            rea = np.pi * 4 * np.real(
                calc_amplitude(
                    kk[0],
                    ampi[3],
                    ampi[2],
                    thet,
                    ampi[1],
                    ampi[0],
                    jset=jselection,
                    jphase=jphas,
                    phase=ph,
                    pola=pola) - calc_amplitude(
                        kk[1],
                        ampi[3],
                        ampi[2],
                        thet,
                        ampi[1],
                        ampi[0],
                        jset=jselection,
                        jphase=jphas,
                        phase=ph,
                        pola=pola))
            ima = np.pi * 4 * np.imag(
                calc_amplitude(
                    kk[0],
                    ampi[3],
                    ampi[2],
                    thet,
                    ampi[1],
                    ampi[0],
                    jset=jselection,
                    jphase=jphas,
                    phase=ph,
                    pola=pola) - calc_amplitude(
                        kk[1],
                        ampi[3],
                        ampi[2],
                        thet,
                        ampi[1],
                        ampi[0],
                        jset=jselection,
                        jphase=jphas,
                        phase=ph,
                        pola=pola))
            tmp.append([rea, ima])
        ampl_of_theta.append(np.array(tmp))

    #for n in [0, int(nbr_angles / 2), -1]:
    #    print(n, ampl_of_theta[n, 0], ampl_of_theta[n, 1], 'I')

    plcols = min(len(indizes), 3)
    plrows = 1 + int(len(indizes) / 3)

    fig = plt.figure(
        num=None, figsize=(15, 12.5), dpi=70, facecolor='w', edgecolor='k')

    plt.title(
        r'$k_1=%2.2f\;MeV, k_2=%2.2f\;MeV$ , $\Delta T:=T^{M_fM_i}_{\lambda^\prime\lambda}(k_1)-T^{M_fM_i}_{\lambda^\prime\lambda}(k_2)$, real (r)BLUE and imaginary (i)RED'
        % (kk[0], kk[1]),
        fontsize=24)

    for k in range(1, 1 + len(indizes)):

        ax1 = fig.add_subplot(plrows - 1, plcols, k)
        ax1.set_xlabel(r'$\theta\;\;[rad]$', fontsize=8)
        ax1.set_ylabel(r'$\Delta T\;\;[fm]$', fontsize=18)

        labre = r'$M_f=%d\;,\;M_i=%d\;,\;\lambda_f=%d\;,\;\lambda_i=%d$' % (
            indizes[k - 1][0], indizes[k - 1][1], indizes[k - 1][2],
            indizes[k - 1][3])
        labim = r''

        ax1.ticklabel_format(axis='y', style='sci', scilimits=(-6, 6))
        ax1.tick_params(axis='y', colors='b')

        # Real
        ax1.plot(
            scattering_angle_range,
            ampl_of_theta[k - 1][:, 0],
            linestyle='--',
            dashes=(5, 2),
            lw=1.5,
            alpha=.85,
            color='b',
            label=labre)
        # Imag
        ax1.plot(
            scattering_angle_range,
            ampl_of_theta[k - 1][:, 1],
            linestyle='--',
            dashes=(5, 8),
            lw=1.5,
            alpha=.95,
            color='r',
            label=labim)

        leg = ax1.legend(loc='best', fontsize=10)
    plt.show()


def plot_crosssection_theta(kk, nbr_angles=40):

    scattering_angle_range = np.linspace(0.0, np.pi, nbr_angles)

    fm2nb = 1e7

    g = {}
    g[0] = np.array((1. / 6.) * (1 + np.cos(scattering_angle_range)**2))
    g[1] = np.array((1. / 4.) * (2 + np.sin(scattering_angle_range)**2))
    g[2] = np.array((1. / 12.) * (13 + np.cos(scattering_angle_range)**2))

    cross_of_theta = np.zeros(nbr_angles)
    for jay in [0, 1, 2]:
        tmp = np.real(
            np.array(
                (polarizability[jay, 'Re'](kk) - 1j * polarizability[jay, 'Im']
                 (kk)) * (polarizability[jay, 'Re']
                          (kk) + 1j * polarizability[jay, 'Im'](kk))))
        cross_of_theta += tmp * g[jay]

    cross_of_theta *= fm2nb

    fig, ax1 = plt.subplots(
        num=None, figsize=(10, 6), dpi=80, facecolor='w', edgecolor='k')
    ax1.set_xlabel(r'$\theta_{cm}\;\;[rad]$', fontsize=8)
    ax1.set_ylabel(r'$(d\sigma/d\Omega)_{cm}\;\;[nb]$', fontsize=8)

    ax1.ticklabel_format(axis='y', style='sci', scilimits=(-6, 6))
    ax1.tick_params(axis='y', colors='b')

    # Real
    ax1.plot(
        scattering_angle_range,
        cross_of_theta,
        linestyle='--',
        dashes=(5, 2),
        lw=1.5,
        alpha=.85,
        color='b',
        label=r'$k_{lab}=%3.1f\,MeV$' % kk)

    #plt.title(r'', fontsize=18)
    leg = ax1.legend(loc='best', fontsize=8)
    plt.show()


# main program section
dbg = 0
jsel = [0, 1, 2]
jphase_exception = []
phas = [1., 1., 1.]  #.076 * np.exp(-1j * np.pi / 15)

MF = -1
MI = -1
pol_phot_out = 1
pol_phot_in = 1
k1 = 40.9746  #MeV
k2 = 10.9746  #MeV

if dbg:

    print('')
    theta = np.pi / 3
    for coss in [[1, 1, 1, 1], [1, 0, 1, 1], [1, -1, 1, 1]]:
        print('(  ', end='')
        for jj in [0, 1, 2]:
            print(
                '%+-4.4f  ' % ((-1)**(coss[3] - coss[1]) * (2 * jj + 1) * w3j(
                    1, jj, 1, -coss[0], coss[0] - coss[1], coss[1]) * w3j(
                        1, 1, jj, coss[2], coss[0] - coss[1] - coss[2],
                        coss[1] - coss[0]) * rotMat(
                            1,
                            (coss[0] - coss[1] - coss[2]), -coss[3], theta)),
                end='')
        print(')')

    jj = 1
    wnk = np.pi / 2.5

    print('d( j=%d  m1=%d  m2=%d )( theta=%2.1f ) = %2.2f' %
          (jj, MF - MI - pol_phot_in, -pol_phot_out, wnk,
           rotMat(jj, MF - MI - pol_phot_in, -pol_phot_out, wnk)))
    print('w3j(1,%d,1,%d,%d,%d) = %4.4f' % (jj, -MF, MF - MI, MI,
                                            w3j(1, jj, 1, -MF, MF - MI, MI)))
    plot_crosssection_theta(kk=k)
    exit()

#plot_pol(2)
#exit()

#print(calc_amplitude(pol_phot_out, 0.0, 'Re', MI, MF))
#plot_ampl(pol_phot_in, pol_phot_out, MI, MF, np.pi / 2)
#exit()

# th extract = Pi/3
hgrie_pol = np.array([
    0.00559711 + 0.00382178 * 1j, 0.00404734 + 0.00550043 * 1j,
    -0.00184007 - 0.00252483 * 1j
]) / (2 * np.pi)

hgrie_pol = []

plot_ampl_of_theta(
    [k1, k2],
    indizes=[[1, 1, 1, 1], [1, 0, 1, 1], [1, -1, 1, 1], [1, 1, 1, -1],
             [1, 0, 1, -1], [1, -1, 1, -1], [1, 1, -1, 1], [1, 0, -1, 1], [
                 1, -1, -1, 1
             ], [1, 1, -1, -1], [1, 0, -1, -1], [1, -1, -1, -1], [0, 0, 1, 1],
             [0, -1, 1, 1], [0, 0, 1, -1], [0, -1, 1, -1], [0, 0, -1, 1],
             [0, -1, -1, 1], [0, 0, -1, -1], [0, -1, -1, -1], [-1, -1, 1, 1], [
                 -1, -1, 1, -1
             ], [-1, -1, -1, 1], [-1, -1, -1, -1]],
    jselection=jsel,
    jphas=jphase_exception,
    ph=phas,
    pola=hgrie_pol)