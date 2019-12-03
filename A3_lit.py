import os, sys
from multiprocessing import Pool

import numpy as np
from scipy.optimize import fmin

from bridge import *
from rrgm_functions import *
from two_particle_functions import *
from three_particle_functions import *
from triton_width_gen import *
from readLITsource_3body import *

os.chdir(v18uixpath)

RHSofBV = {}
RHSofmJ = {}

if 'construe_fresh_helion' in cal:

    os.system('cp INLUCN_avui INLUCN')
    os.system(BINBDGpath + 'LUDW_CN.exe')
    os.system('cp INLU_avui INLU')
    os.system(BINBDGpath + 'DRLUD.exe')
    os.system('cp INOB_avui INOB')
    os.system(BINBDGpath + 'KOBER.exe')
    os.system(BINBDGpath + 'DROBER.exe')
    os.system('cp INQUA_N_V18 INQUA_N')
    os.system(BINBDGpath + 'QUAFL_N.exe')
    os.system('cp INQUA_N_UIX INQUA_N')
    os.system(BINBDGpath + 'DRQUA_AK_N.exe')
    os.system('cp INEN_avui INEN')
    os.system(BINBDGpath + 'DR2END_AK.exe')

    EBDG = get_h_ev()[0]
    np.savetxt('E0.dat', np.array([EBDG]), fmt='%12.4f')

    os.system('cp OUTPUT end_out_b && cp INEN inen_b')
    rrgm_functions.parse_ev_coeffs(infil='end_out_b')

if 'rhs' in cal:

    for file in os.listdir(litpath3He):
        if fnmatch.fnmatch(file, 'endlit*'):
            if 'dbg' in cal:
                print('removing old <*en*lit*> files.')
            os.system('rm ' + litpath3He + '/*en*lit*')
            break

    for streukanal in streukas:

        os.chdir(v18uixpath)

        Jstreu = float(streukanal.split('^')[0])
        mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0,
                                                      Jstreu)

        BUECO = [cof.strip() for cof in open('COEFF')]
        EBDG = get_h_ev(ifi='end_out_b')[0]

        print('(iv)    LS-scheme: B(2,%s) = %4.4f MeV [' % (boundstatekanal,
                                                            EBDG),
              get_h_ev(n=4, ifi='end_out_b'), ']')
        print('        dim(B_0)   = %d' % len(BUECO))

        os.chdir(litpath3He)
        # assume: EINE Zerlegung pro L und S cfg. => anz int.ws < x

        lfrags = []
        sfrags = []
        anzBSzerl = 0
        for lcfg in range(len(rechtekanaele[boundstatekanal])):
            sfrags = sfrags + rechtekanaele[boundstatekanal][lcfg][1]
            for scfg in rechtekanaele[boundstatekanal][lcfg][1]:
                anzBSzerl += 1
                lfrags = lfrags + [rechtekanaele[boundstatekanal][lcfg][0]]
        for lcfg in range(len(streukanaele3He[streukanal])):
            sfrags = sfrags + streukanaele3He[streukanal][lcfg][1]
            for scfg in streukanaele3He[streukanal][lcfg][1]:
                lfrags = lfrags + [streukanaele3He[streukanal][lcfg][0]]

        if 'dbg' in cal:
            print(lfrags, sfrags)

        if 'lit_lu-ob-qua' in cal:

            lit_3inlu(mul=multipolarity, frag=lfrags)
            os.system(BINLITpath + 'luise.exe > dump')
            lit_3inob(fr=sfrags)
            os.system(BINLITpath + 'obem.exe > dump')

            intw = threebodywidths(sfrags[anzBSzerl:], nPERz=basisdimLITint)
            np.savetxt('intw3heLIT.dat', np.array(intw), fmt='%12.4f')
            relw = threebodywidths(sfrags[anzBSzerl:], nPERz=basisdimLITrel)
            np.savetxt('relw3heLIT.dat', np.array(relw), fmt='%12.4f')
            if 'dbg' in cal:
                print('LIT-basis widths (internal):\n', intw)
                print('LIT-basis widths (relative):\n', relw)

            lit_3inqua(anzo=11, bnd=v18uixpath + '/INQUA_N')
            lit_3inqua(intwi=intw, relwi=relw, withhead=False)

            os.system(BINLITpath + 'qual.exe')

        leftpar = 1 if streukanal[-1] == '-' else 2

        def cal_rhs_compo(para):
            pid = os.getpid()
            print('process id:', pid)
            os.system('mkdir tmp.' + str(pid))
            os.chdir('tmp.' + str(pid))
            os.system('cp ../QUAOUT .')

            lit_3inen(
                MREG='  1  0  0  0  0  0  0  0  0  1  1',
                #                   (shifted) QBV                     nr.rw
                KSTREU=[int(para[1] + len(BUECO)), para[2]],
                JWSL=Jstreu,
                JWSLM=para[0][1],
                MULM2=para[0][0],
                NPARL=leftpar,
                JWSR=J0,
                NPARR=2,
                EB=EBDG,
                BUECO=BUECO,
                NZE=anz_phot_e,
                EK0=phot_e_0,
                EKDIFF=phot_e_d,
                bnd=v18uixpath + '/INEN_avui')

            print(para)
            os.system(BINLITpath + 'enemb.exe')
            os.system('cp OUTPUT ../endlit%d-%d_J%f_mJ%f-mL%d.log' %
                      (para[1], para[2], Jstreu, para[0][1], para[0][0]))
            os.system('cp INEN ../inenlit%d-%d_J%f_mJ%f-mL%d.log' %
                      (para[1], para[2], Jstreu, para[0][1], para[0][0]))
            os.chdir(litpath3He)
            os.system('rm -rf tmp.' + str(pid))

        parameter_set = []
        litdim = sum([len(sc[1]) for sc in streukanaele3He[streukanal]
                      ]) * basisdimLITint * basisdimLITrel

        print(
            'N = %d basis vectors can be combined with %d relative widths\n => |LIT basis| = %d'
            % (sum([len(sc[1])
                    for sc in streukanaele3He[streukanal]]) * basisdimLITint,
               basisdimLITrel, litdim))

        # ALL [ BV nr , relw ] tuples
        bv_rel_combos = np.array(
            np.meshgrid(
                np.arange(1, 1 + sum(
                    [len(sc[1])
                     for sc in streukanaele3He[streukanal]]) * basisdimLITint),
                np.arange(1, 1 + basisdimLITrel))).T.reshape(-1, 2)

        # random tuple subset of dinemsion LD
        idx = np.random.choice(len(bv_rel_combos), LD)
        litbas = np.unique(bv_rel_combos[idx], axis=0)
        np.savetxt('LITbas.dat', litbas, fmt='%d')

        for mM in mLmJl:
            bvnstreu = 0
            for bv in litbas:
                bvnstreu += 1
                parameter_set.append([mM, bv[0], bv[1]])

        def pool_meister():
            p = Pool(anzproc)
            p.map(cal_rhs_compo, parameter_set)

        pool_meister()
        exit()

if 'lhs' in cal:

    print('(ii)    calculating norm/ham in scattering-channel basis')
    os.chdir(v18uixpath)

    for streukanal in streukas:
        lfrags = []
        sfrags = []
        anzBSzerl = 0
        for lcfg in range(len(streukanaele3He[streukanal])):
            sfrags = sfrags + streukanaele3He[streukanal][lcfg][1]
            for scfg in streukanaele3He[streukanal][lcfg][1]:
                lfrags = lfrags + [streukanaele3He[streukanal][lcfg][0]]

        print(lfrags, sfrags)
        Jstreu = float(streukanal.split('^')[0])

        n3_inlu(21, fn='INLU', fr=lfrags)
        os.system(BINBDGpath + 'DRLUD.exe')
        n3_inlu(21, fn='INLUCN', fr=lfrags)
        os.system(BINBDGpath + 'LUDW_CN.exe')

        n3_inob(sfrags, 20, fn='INOB')
        os.system(BINBDGpath + 'KOBER.exe')
        os.system(BINBDGpath + 'DROBER.exe')

        intwLIT = np.loadtxt(litpath3He + 'intw3heLIT.dat').reshape(
            (len(sfrags), -1))
        anzLITbv = sum([len(frgm) for frgm in intwLIT])
        relwLIT = np.loadtxt(litpath3He + 'relw3heLIT.dat').reshape(
            (len(sfrags), -1))
        he3inqua(
            intwi=intwLIT,
            relwi=relwLIT,
            potf='/home/kirscher/kette_repo/sim_par/potentials/NN_pheno/AV18')
        os.system(BINBDGpath + 'QUAFL_N.exe')
        repl_line(
            'INQUA_N', 1,
            '/home/kirscher/kette_repo/sim_par/potentials/NNN_pheno/urbana9_AK_neu\n'
        )
        os.system(BINBDGpath + 'DRQUA_AK_N.exe')

        costr = ''
        for nn in range(1, 30):
            cf = 1.0 if nn < 28 else 0.0
            costr += '%12.7f' % cf if (nn % 7 != 0) else '%12.7f\n' % cf

        litbas = np.loadtxt(litpath3He + 'LITbas.dat').astype(int)
        n3_inen_rhs(
            litbas,
            Jstreu,
            costr,
            np.ones(len(relwLIT[0])),
            fn='INEN',
            pari=0,
            nzop=31,
            tni=11)

        os.system(BINBDGpath + 'DR2END_AK.exe')

        os.system('cp INEN inen-lit-%s' % streukanal)
        os.system('cp %s/MATOUT %s/norm-ham-litME-%s' %
                  (v18uixpath, v18uixpath, streukanal))

if 'plt' in cal:

    litbas = np.loadtxt(litpath3He + 'LITbas.dat').astype(int)[:3]

    for streukanal in streukas:
        # read uncoupled source ME's
        os.chdir(litpath3He)
        RHSofBV[streukanal], photEn = read_uncoupled_source(
            streukanal, basisSET=litbas)
        # couple incoming state with photon multipole to Jlit
        exit()
        RHSofmJ[streukanal] = couple_source(
            streukanal, RHSofBV[streukanal], basisSET=litbas)

    fig = plt.figure()
    fig.subplots_adjust(hspace=0.4, wspace=0.4)
    for i in range(len(streukas)):
        Jstreu = float(streukas[i].split('^')[0])
        ax1 = fig.add_subplot(len(streukas), 2, 2 * i + 1)
        #ax1.set_title(r'$J^\pi=%d^%s$' % (Jstreu, streukas[i][1]))
        ax1.set_xlabel('photon momentum [MeV]')
        #ax1.set_title(r'$J^\pi=%d^%s$' % (Jstreu, streukas[i][1]))
        mLmJl, mLrange, mJlrange = non_zero_couplings(multipolarity, J0,
                                                      Jstreu)
        mM = mLmJl[0]
        [
            ax1.plot(photEn,
                     RHSofBV[streukas[i]][('%d-%d' % (bv[0], bv[1]), '%d' %
                                           (2 * Jstreu), '%d' % (2 * mM[1]),
                                           '%d' % (2 * multipolarity),
                                           '%d' % (2 * mM[0]))])
            for bv in litbas
        ]

        #ax.text(0.5, 0.5, str((2, 3, i)),
        #       fontsize=18, ha='center')
        ax2 = fig.add_subplot(len(streukas), 2, 2 * i + 2)
        ax2.set_xlabel('photon momentum [MeV]')
        ax2.set_ylabel(r'$\left\langle\,Jm\,\vert\,Jm\,\right\rangle$ [-]')
        ax2.set_title(r'$J$-coupled RHS')

        [
            ax2.plot(photEn,
                     RHSofmJ[streukas[i]][('%d-%d' % (bv[0], bv[1]), '%d' %
                                           (2 * Jstreu), '%d' % (2 * mM[1]),
                                           '%d' % (2 * multipolarity))])
            for bv in litbas
        ]

    #plt.show()
    fig.savefig(v18uixpath + '/LITrhs.pdf')