import os, sys

import numpy as np
from scipy.optimize import fmin

from bridge import *
from rrgm_functions import *
from two_particle_functions import *
from lit_plots import *
from readLITsource import *

from plot_basis import main

costr = ''
for nn in range(14):
    costr += '%12.7f' % 1.0 if (nn != 6) else '%12.7f\n' % 1.0

os.chdir(av18path)
RHSofBV = {}
RHSofmJ = {}

if 'construe_fresh_deuteron' in cal:

    h2_inlu(nfrag=nzf0)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    h2_inob(nfrag=nzf0)
    os.system(BINBDGpath + 'KOBER.exe')
    h2_inqua(rw0, potnn)
    os.system(BINBDGpath + 'QUAFL_N.exe')

    h2_inen_bs(
        relw=rw0,
        costr=costr,
        j=int(boundstatekanal[-1]),
        ch=rechtekanaele[boundstatekanal],
        nfrag=nzf0)

    os.system(BINBDGpath + 'DR2END_AK.exe')
    exit()
    purge_basis(
        max_coeff=10000,
        min_coeff=150,
        nbr_cycles=40,
        max_diff=0.05,
        dr2executable=BINBDGpath + 'DR2END_AK.exe',
        dbg=True if ('dbg' in cal) else False)

    clean_wrels, clean_inen_indices = purged_width_set(infil='INEN', w0=rw0)
    np.savetxt('w0.dat', np.array(clean_wrels), fmt='%12.4f')

    nzf0p = int(np.ceil(len(clean_wrels) / 20.0))

    if (nzf0 != nzf0p):
        print('NZF0 (after purge) != NZF0 (before purge)', nzf0, nzf0p)
        h2_inlu(nfrag=nzf0p)
        os.system(BINBDGpath + 'LUDW_CN.exe')
        h2_inob(nfrag=nzf0p)
        os.system(BINBDGpath + 'KOBER.exe')

    h2_inqua(clean_wrels, potnn)
    for ind_set in range(len(clean_inen_indices)):
        istr = ''
        for i in clean_inen_indices[ind_set]:
            istr += '%3s' % str(i)
        istr += '\n'
        repl_line('INEN', int(6 + 2 * ind_set), istr)

    os.system(BINBDGpath + 'QUAFL_N.exe')
    os.system(BINBDGpath + 'DR2END_AK.exe')
    EBDG = get_h_ev()[0]
    np.savetxt('E0.dat', np.array([EBDG]), fmt='%12.4f')

    os.system('cp OUTPUT end_out_b && cp INEN inen_b')
    rrgm_functions.parse_ev_coeffs(infil='end_out_b')

if 'construe_fresh_LIT_basis' in cal:

    clean_wrels = [float(line) for line in open('w0.dat')]
    h2_inqua(clean_wrels, potnn)
    nzf0p = int(np.ceil(len(clean_wrels) / 20.0))

    wLIT = sparsifyOnlyOne(winiLIT, clean_wrels, min_spacing_to_LITWs)
    np.savetxt('wLIT.dat', np.array(wLIT), fmt='%12.4f')

    nzfLIT = int(np.ceil(len(wLIT) / 20.0))
    nzfTOT = nzf0p + nzfLIT

    # consistency check if B0 is the same as in the smaller space
    h2_inlu(nfrag=nzfTOT)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    h2_inob(nfrag=nzfTOT)
    os.system(BINBDGpath + 'KOBER.exe')
    h2_inqua(wLIT, potnn, withhead=False)
    os.system(BINBDGpath + 'QUAFL_N.exe')
    os.system('cp inen_b INEN')
    os.system(BINBDGpath + 'DR2END_AK.exe')

    EBDG = get_h_ev()[0]
    E0small = [float(line) for line in open('E0.dat')][0]
    print('%8.4f == %8.4f ?' % (EBDG, E0small))

    main()

for file in os.listdir(litpathD):
    if fnmatch.fnmatch(file, 'endlit*'):
        if 'dbg' in cal:
            print('removing old <*en*lit*> files.')
        os.system('rm ' + litpathD + '/*en*lit*')
        break

for streukanal in streukas:

    os.chdir(av18path)

    Jstreu = float(streukanal.split('^')[0])
    mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0, Jstreu)

    wLIT = np.loadtxt('wLIT.dat')
    clean_wrels = np.loadtxt('w0.dat')
    nzfLIT = int(np.ceil(len(wLIT) / 20.0))
    nzf0p = int(np.ceil(len(clean_wrels) / 20.0))
    nzfTOT = nzf0p + nzfLIT

    BUECO = [cof.strip() for cof in open('COEFF')]
    BSRWIDX = get_bsv_rw_idx(chs=2, inen='inen_b')
    EBDG = get_h_ev(ifi='end_out_b')[0]

    print(
        '(iv)    LS-scheme: B(2,%s) = %4.4f MeV [' % (boundstatekanal, EBDG),
        get_h_ev(n=4),
        ']')
    print('        dim(B_0)   = %d -> %d' % (len(rw0), len(clean_wrels)))
    print('        dim(B_LIT) = %d -> %d' % (len(winiLIT), len(wLIT)))
    print('        dim(B)     = %d' % (len(wLIT) + len(clean_wrels)))

    os.chdir(litpathD)

    lit_inlu(mul=multipolarity, nfrag=nzfTOT)
    os.system(BINLITpath + 'luise.exe > dump')
    lit_inob(nfrag=nzfTOT)
    os.system(BINLITpath + 'obem.exe > dump')

    lit_inqua(anzo=11, relw=clean_wrels)
    lit_inqua(relw=wLIT, withhead=False, LREG='')

    os.system(BINLITpath + 'qual.exe')

    leftpar = 1 if streukanal[-1] == '-' else 2

    for mM in mLmJl:
        for subchannel in range(len(streukanaeleD[streukanal])):
            for streukanalweite in range(1, len(wLIT) + 1):
                lit_inen(
                    MREG='  1  0  0  0  0  0  0  0  0  1  1',
                    #                   (shifted) QBV                     nr.rw
                    KSTREU=[
                        int(streukanaeleD[streukanal][subchannel][-1] +
                            8 * nzf0p), streukanalweite
                    ],
                    JWSL=Jstreu,
                    JWSLM=mM[1],
                    MULM2=mM[0],
                    NPARL=leftpar,
                    KBND=[[2, BSRWIDX[0]], [6, BSRWIDX[1]]],
                    JWSR=1,
                    NPARR=2,
                    EB=EBDG,
                    BUECO=BUECO,
                    NZE=anz_phot_e,
                    EK0=phot_e_0,
                    EKDIFF=phot_e_d)

                os.system(BINLITpath + 'enemb.exe')
                os.system('cp OUTPUT endlit%d_J%d_mJ%d-mL%d.log' %
                          (int(streukanalweite + subchannel * len(wLIT)),
                           Jstreu, mM[1], mM[0]))
                os.system('cp INEN inenlit%d_J%d_mJ%d-mL%d.log' %
                          (int(streukanalweite + subchannel * len(wLIT)),
                           Jstreu, mM[1], mM[0]))

    if 'dbg' in cal:
        print('(iiib)  calculated S_bv^(Jlit,m) for mL,m in:', mLmJl)

    print('(ii)    calculating norm/ham in scattering-channel basis')
    os.chdir(av18path)

    h2_inlu(nfrag=nzfLIT)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    h2_inob(nfrag=nzfLIT)
    os.system(BINBDGpath + 'KOBER.exe')
    h2_inqua(wLIT, potnn)
    os.system(BINBDGpath + 'QUAFL_N.exe')

    h2_inen_bs(
        relw=wLIT,
        costr=costr,
        j=Jstreu,
        ch=streukanaeleD[streukanal],
        nfrag=nzfLIT)

    os.system(BINBDGpath + 'DR2END_AK.exe')

    os.system('cp INEN inen-lit-%s' % streukanal)
    os.system('cp %s/MATOUT %s/norm-ham-litME-%s' % (av18path, av18path,
                                                     streukanal))

    # read uncoupled source ME's
    os.chdir(litpathD)
    RHSofBV[streukanal], photEn = read_uncoupled_source(
        streukanal, basisSET=wLIT)
    # couple incoming state with photon multipole to Jlit
    RHSofmJ[streukanal] = couple_source(
        streukanal, RHSofBV[streukanal], basisSET=wLIT)

fig = plt.figure()
fig.subplots_adjust(hspace=0.4, wspace=0.4)
for i in range(len(streukas)):
    ax1 = fig.add_subplot(len(streukas), 2, 2 * i + 1)
    ax1.set_title(r'$J^\pi=%d^%s$' % (int(streukas[i][0]), streukas[i][1]))
    ax1.set_xlabel('photon momentum [MeV]')
    ax1.set_title(r'$J^\pi=%d^%s$' % (int(streukas[i][0]), streukas[i][1]))
    mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0,
                                                  int(streukas[i][0]))
    mM = mLmJl[0]
    [
        ax1.plot(photEn,
                 RHSofBV[streukas[i]][('%d' % (streukanalweite),
                                       '%d' % (2 * int(streukas[i][0])), '%d' %
                                       (2 * mM[1]), '%d' % (2 * multipolarity),
                                       '%d' % (2 * mM[0]))])
        for streukanalweite in range(1,
                                     len(wLIT) + 1)
    ]

    #ax.text(0.5, 0.5, str((2, 3, i)),
    #       fontsize=18, ha='center')
    ax2 = fig.add_subplot(len(streukas), 2, 2 * i + 2)
    ax2.set_xlabel('photon momentum [MeV]')
    ax2.set_ylabel(r'$\left\langle\,Jm\,\vert\,Jm\,\right\rangle$ [-]')
    ax2.set_title(r'$J$-coupled RHS')

    [
        ax2.plot(photEn, RHSofmJ[streukas[i]][('%d' % (streukanalweite), '%d' %
                                               (2 * int(streukas[i][0])),
                                               '%d' % (2 * mM[1]),
                                               '%d' % (2 * multipolarity))])
        for streukanalweite in range(1,
                                     len(wLIT) + 1)
    ]

#plt.show()
fig.savefig(av18path + '/LITrhs.pdf')