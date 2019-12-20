from bridge import *
from three_particle_functions import *
from triton_width_gen import *
import operator

os.chdir(v18uixpath)

for streukanal in streukas:

    lfrags = []
    sfrags = []

    for lcfg in range(len(streukanaele3He[streukanal])):
        sfrags = sfrags + streukanaele3He[streukanal][lcfg][1]
        for scfg in streukanaele3He[streukanal][lcfg][1]:
            lfrags = lfrags + [streukanaele3He[streukanal][lcfg][0]]

    print('\n3-body LIT-basis components (bare):%5d :' % len(lfrags))

    # read width sets of the ground-state (3-helium) basis
    he_iw, he_rw, frgs = retrieve_he3_widths(v18uixpath + 'INQUA_N_UIX')

    # generate widths for the two coordinates for every LIT-basis spin-/angular-momentum comfiguration
    seedwidths1 = [
        10.9702, 7.2388, 0.7274, 0.1739, 0.03548, 2.8716, 0.9217, 0.2435, 3.5,
        4, 4.5, 5.5
    ]
    seedwidths2 = seedwidths1  #np.sort(np.unique(np.random.choice(w12, 15)))[::-1]

    lfrags2 = []
    sfrags2 = []
    min_single_width_dist = 0.15
    min_eucl_pair_dist = .2
    lit_iw = []
    lit_rw = []

    widthcluster = 15

    for frg in range(len(lfrags)):

        rho1widths = perturbseedwidths(seedwidths1, widthcluster).tolist()
        rho2widths = perturbseedwidths(seedwidths2, widthcluster).tolist()

        iw = sparsify(rho1widths, min_single_width_dist)
        rw = sparsify(rho2widths, min_single_width_dist)

        lit_iw += np.array(np.array_split(iw, int(np.ceil(
            len(iw) / 12.)))).tolist()

        lit_rw += int(np.ceil(len(iw) / 12.)) * rw.reshape((1, -1)).tolist()

        sfrags2 += int(np.ceil(len(iw) / 12.)) * [sfrags[frg]]
        lfrags2 += int(np.ceil(len(iw) / 12.)) * [lfrags[frg]]

    with open(litpath3He + 'frags_LIT.dat', 'wb') as f:
        np.savetxt(
            f,
            np.column_stack([sfrags2, lfrags2]),
            fmt='%s',
            delimiter=' ',
            newline=os.linesep)
    f.close()

    with open(litpath3He + 'intw3heLIT.dat', 'wb') as f:
        for ws in lit_iw:
            np.savetxt(f, [ws], fmt='%12.4f', delimiter=' ; ')
    f.close()
    with open(litpath3He + 'relw3heLIT.dat', 'wb') as f:
        for ws in lit_rw:
            np.savetxt(f, [ws], fmt='%12.4f', delimiter=' ; ')
    f.close()

    bas = []
    for nz in range(len(lit_iw)):
        allindices = [[w1, w2] for w1 in range(len(lit_iw[nz]))
                      for w2 in range(len(lit_rw[nz]))]
        e0 = np.random.choice(range(len(allindices)))
        set0 = [allindices[e0]]
        allindices.remove(allindices[e0])
        while len(allindices) > 0:
            e0 = np.random.choice(range(len(allindices)))

            reject = False
            for ei in range(len(set0)):

                w1 = [
                    lit_iw[nz][allindices[e0][0]],
                    lit_rw[nz][allindices[e0][1]]
                ]
                w2 = [lit_iw[nz][set0[ei][0]], lit_rw[nz][set0[ei][1]]]

                if np.linalg.norm(
                        np.array(w1) - np.array(w2)) < min_eucl_pair_dist:
                    reject = True
                    break

            if reject:
                allindices.remove(allindices[e0])
                continue
            else:
                set0.append(allindices[e0])
                allindices.remove(allindices[e0])

        offset_int = sum([len(iws) for iws in lit_iw[:nz]])
        bas += [[el[0] + 1 + offset_int, el[1] + 1] for el in set0]

    basvs = {}
    for bv in bas:
        try:
            basvs[bv[0]].append(bv[1])
        except:
            basvs[bv[0]] = [bv[1]]
    sbas = sorted(basvs.items(), key=operator.itemgetter(0))

    with open(v18uixpath + 'LITbas_full.dat', 'wb') as f:
        np.savetxt(f, [[jj[0], kk] for jj in sbas for kk in jj[1]], fmt='%d')
    f.close()

    n3_inlu(21, fn='INLU', fr=lfrags2)
    os.system(BINBDGpath + 'DRLUD.exe')
    n3_inlu(21, fn='INLUCN', fr=lfrags2)
    os.system(BINBDGpath + 'LUDW_CN.exe')
    n3_inob(sfrags2, 20, fn='INOB')
    os.system(BINBDGpath + 'KOBER.exe')
    os.system(BINBDGpath + 'DROBER.exe')
    he3inqua(
        intwi=lit_iw,
        relwi=lit_rw,
        potf='/home/kirscher/kette_repo/sim_par/potentials/NN_pheno/AV18')
    os.system('time ' + BINBDGpath + 'QUAFL_N.exe')
    repl_line(
        'INQUA_N', 1,
        '/home/kirscher/kette_repo/sim_par/potentials/NNN_pheno/urbana9_AK_neu\n'
    )
    os.system('time ' + BINBDGpath + 'DRQUA_AK_N.exe')

    n3_inen_bdg(sbas, 1. / 2., costr, fn='INEN', pari=0, nzop=31, tni=11)

    os.system(BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')