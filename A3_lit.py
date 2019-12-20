import os, sys
import multiprocessing
import subprocess
import shlex
from multiprocessing.pool import ThreadPool

import numpy as np
from scipy.optimize import fmin

from bridge import *
from rrgm_functions import *
from two_particle_functions import *
from three_particle_functions import *
from triton_width_gen import *
from readLITsource_3body import *
from BasisVisualization import visbas

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

    print('helion ground state calculated!')

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

        if 'rhs_lu-ob-qua' in cal:

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

            # read width sets of the ground-state (3-helium) basis
            he_iw, he_rw, frgs = retrieve_he3_widths(
                v18uixpath + 'INQUA_N_UIX')

            # read widths and frags of the LIT basis as determined via
            # v18uix_LITbasis.py
            fragfile = [ln for ln in open(litpath3He + 'frags_LIT.dat')]
            lfrags2 = [fr.split(' ')[1].strip() for fr in fragfile]
            sfrags2 = [fr.split(' ')[0] for fr in fragfile]

            intwLIT = [
                np.array(ln.split(';')).astype(float).tolist()
                for ln in open(litpath3He + 'intw3heLIT.dat')
            ]

            relwLIT = [
                np.array(ln.split(';')).astype(float).tolist()
                for ln in open(litpath3He + 'relw3heLIT.dat')
            ]

            if 'dbg' in cal:
                print(
                    '\n3He components (full) + LIT-basis components (bare):\n',
                    len(lfrags))
                print(sfrags)
                print('\nLIT-basis components (full):\n', len(lfrags2))
                print(sfrags2)

            lit_3inqua(
                anzo=11,
                LREG='  1  0  0  0  0  0  0  0  0  1  1',
                bnd=v18uixpath + '/INQUA_N_UIX',
                outfile='INQUA_3he')
            lit_3inqua(
                intwi=intwLIT,
                relwi=relwLIT,
                withhead=False,
                outfile='INQUA_LIT')

            os.system('cat INQUA_3he INQUA_LIT > INQUA')
            os.system('cp INQUA INQUA_3he+LIT')

            lit_3inlu(mul=multipolarity, frag=lfrags[:anzBSzerl] + lfrags2)
            os.system(BINLITpath + 'luise.exe > dump')
            lit_3inob(fr=sfrags[:anzBSzerl] + sfrags2)
            os.system(BINLITpath + 'obem.exe > dump')

            with open('frags_3He+LIT.dat', 'wb') as f:
                np.savetxt(
                    f,
                    np.column_stack([
                        ['#BV=%d' % anzBSzerl] + sfrags[:anzBSzerl] + sfrags2,
                        ['#BV=%d' % anzBSzerl] + lfrags[:anzBSzerl] + lfrags2
                    ]),
                    fmt='%s',
                    delimiter=' ',
                    newline=os.linesep)
            f.close()

            os.system(BINLITpath + 'qual.exe')

        leftpar = 1 if streukanal[-1] == '-' else 2

        def cal_rhs_compo(para, procnbr):

            inenf = 'inenlit%d-%d_J%f_mJ%f-mL%d.log' % (para[1], para[2],
                                                        Jstreu, para[0][1],
                                                        para[0][0])
            outfs = 'endlit%d-%d_J%f_mJ%f-mL%d.log FORMFA_%s' % (para[1],
                                                                 para[2],
                                                                 Jstreu,
                                                                 para[0][1],
                                                                 para[0][0],
                                                                 str(procnbr))

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
                bnd=v18uixpath + '/INEN_avui',
                outfile=inenf)

            cmd = BINLITpath + 'enemb.exe %s %s' % (inenf, outfs)
            p = subprocess.Popen(
                shlex.split(cmd),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            out, err = p.communicate()
            return (out, err)

        parameter_set = []

        if 'new_rnd_bvset' in cal:
            intwLIT = [
                np.array(ln.split(';')).astype(float)
                for ln in open(litpath3He + 'intw3heLIT.dat')
            ]

            basisdimLIT = sum([len(ws) for ws in intwLIT])

            anzLITbv = sum([len(frgm) for frgm in intwLIT])

            relwLIT = [
                np.array(ln.split(';')).astype(float)
                for ln in open(litpath3He + 'relw3heLIT.dat')
            ]

            litdim = 0
            anzbv = 0
            anzrel = 0
            bv_rel_combos = np.empty([1, 2])
            for z in range(len(intwLIT)):
                anzbv += len(intwLIT[z])
                anzrel += len(relwLIT[z])
                litdim += len(intwLIT[z]) * len(relwLIT[z])
                # ALL [ BV nr , relw ] tuples
                bv_rel_combos = np.concatenate(
                    (bv_rel_combos,
                     np.array(
                         np.meshgrid(
                             np.arange(
                                 1 + sum([len(iss) for iss in intwLIT[:z]]),
                                 1 + sum([len(iss) for iss in intwLIT[:z]]) +
                                 len(intwLIT[z])),
                             np.arange(1, 1 + len(relwLIT[z])))).T.reshape(
                                 -1, 2)),
                    axis=0).astype(int)

            bv_rel_combos = bv_rel_combos[1:]

            print(
                'N = %d basis vectors can be combined with %d relative widths\n => |LIT basis| = %d'
                % (anzbv, anzrel, litdim))
            # random tuple subset of dinemsion LD
            litbas = []
            while (len(litbas) < LD):
                idx = np.random.choice(len(bv_rel_combos), 4)
                litbas += [
                    tupel
                    for tupel in list(np.unique(bv_rel_combos[idx], axis=0))
                ]
                litbas = list(np.unique(litbas, axis=0))

            print('(randomly) selected BV-rel tuples:\n',
                  [list(jj) for jj in litbas])

            with open(v18uixpath + 'LITbas_rnd.dat', 'wb') as f:
                np.savetxt(f, litbas, fmt='%d')
            f.close()

        else:
            litbas = np.loadtxt(v18uixpath + 'LITbas_red.dat').astype(int)
            print('[...] reading BV-rw tupel from <LITbas_red.dat>:\n',
                  [list(jj) for jj in litbas])
            visbas(
                basispath=v18uixpath, widthpath=litpath3He, exepath=BINBDGpath)
            os.chdir(litpath3He)
            print('basis width grid plotted in <WidthXY.pdf>')

        for mM in mLmJl:
            bvnstreu = 0
            for bv in litbas:
                bvnstreu += 1
                parameter_set.append([mM, bv[0], bv[1]])

        pool = ThreadPool(anzproc)

        results = []

        for procnbr in range(len(parameter_set)):
            pars = parameter_set[procnbr]
            results.append(pool.apply_async(cal_rhs_compo, (
                pars,
                procnbr,
            )))

        pool.close()
        pool.join()

if 'lhs' in cal:

    print('(ii)    calculating norm/ham in scattering-channel basis')
    os.chdir(v18uixpath)

    for streukanal in streukas:

        fragfile = [ln for ln in open(litpath3He + 'frags_3He+LIT.dat')]
        anzF3he = int(fragfile[0].strip()[-1])
        lfrags = [fr.split(' ')[1].strip() for fr in fragfile[1 + anzF3he:]]
        sfrags = [fr.split(' ')[0] for fr in fragfile[1 + anzF3he:]]

        intwLIT = [
            np.array(ln.split(';')).astype(float)
            for ln in open(litpath3He + 'intw3heLIT.dat')
        ]

        anzLITbv = sum([len(frgm) for frgm in intwLIT])

        if (len([len(frgm) for frgm in intwLIT]) !=
            (len(fragfile) - 1 - anzF3he)):
            print('LIT-basis fragments inconcistent!',
                  len([len(frgm) for frgm in intwLIT]),
                  (len(fragfile) - 1 - anzF3he))
            exit()

        relwLIT = [
            np.array(ln.split(';')).astype(float)
            for ln in open(litpath3He + 'relw3heLIT.dat')
        ]

        if 'dbg' in cal:
            print(lfrags, sfrags)
            #print(intwLIT, relwLIT)

        he_iw, he_rw, frags = retrieve_he3_widths(litpath3He + 'INQUA_LIT')

        Jstreu = float(streukanal.split('^')[0])

        if 'lhs_lu-ob-qua' in cal:

            n3_inlu(21, fn='INLU', fr=lfrags)
            os.system(BINBDGpath + 'DRLUD.exe')
            n3_inlu(21, fn='INLUCN', fr=lfrags)
            os.system(BINBDGpath + 'LUDW_CN.exe')

            n3_inob(sfrags, 20, fn='INOB')
            os.system(BINBDGpath + 'KOBER.exe')
            os.system(BINBDGpath + 'DROBER.exe')

            he3inqua(
                intwi=intwLIT,
                relwi=relwLIT,
                potf=
                '/home/kirscher/kette_repo/sim_par/potentials/NN_pheno/AV18')
            os.system(BINBDGpath + 'QUAFL_N.exe')
            repl_line(
                'INQUA_N', 1,
                '/home/kirscher/kette_repo/sim_par/potentials/NNN_pheno/urbana9_AK_neu\n'
            )
            os.system(BINBDGpath + 'DRQUA_AK_N.exe')

        litbas = np.loadtxt(v18uixpath + 'LITbas_red.dat').astype(int)

        n3_inen_rhs(
            litbas,
            Jstreu,
            costr,
            np.ones(len(relwLIT[0])),
            fn='INEN',
            pari=0,
            nzop=31,
            tni=11)

        os.system(
            BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')

        os.system('cp INEN inen-lit-%s' % streukanal)
        os.system('cp %s/MATOUT %s/norm-ham-litME-%s' %
                  (v18uixpath, v18uixpath, streukanal))

if 'plt' in cal:

    litbas = np.loadtxt(v18uixpath + 'LITbas_red.dat').astype(int)

    for streukanal in streukas:
        # read uncoupled source ME's
        os.chdir(litpath3He)
        RHSofBV[streukanal], photEn = read_uncoupled_source(
            streukanal, basisSET=litbas)
        # couple incoming state with photon multipole to Jlit
        RHSofmJ[streukanal] = couple_source(
            streukanal, RHSofBV[streukanal], basisSET=litbas)

    fig = plt.figure(figsize=(10, 8), dpi=95)

    fig.subplots_adjust(hspace=1.4, wspace=0.4)
    for i in range(len(streukas)):
        Jstreu = float(streukanal.split('^')[0])
        ax1 = fig.add_subplot(len(streukas), 2, 2 * i + 1)
        ax1.set_title(r'$J^\pi=%d^%s$' % (Jstreu, streukas[i][-1]))
        ax1.set_xlabel('photon momentum [MeV]')
        #ax1.set_title(r'$J^\pi=%d^%s$' % (Jstreu, streukas[i][-1]))
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
            ax2.plot(
                photEn,
                RHSofmJ[streukas[i]][('%d-%d' % (bv[0], bv[1]),
                                      '%d' % (2 * Jstreu), '%d' % (2 * mM[1]),
                                      '%d' % (2 * multipolarity))],
                label=r'$BV_{%d}W^{%d}$' % (int(bv[0]), int(bv[1])))
            for bv in litbas
        ]

        box1 = ax1.get_position()
        ax1.set_position([box1.x0, box1.y0, box1.width, box1.height * 0.5])
        box2 = ax2.get_position()
        ax2.set_position([box2.x0, box2.y0, box2.width, box2.height * 0.5])

        leg = ax2.legend(
            loc='upper center',
            ncol=8,
            numpoints=1,
            fontsize=10,
            bbox_to_anchor=(-0.25, 2.3))

    #plt.show()
    fig.savefig(v18uixpath + 'LITrhs.pdf')