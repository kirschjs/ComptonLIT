from bridge import *
from three_particle_functions import *
from BasisVisualization import visbas
import operator

os.chdir(v18uixpath)
print('(working dir) %s' % v18uixpath)

bas = np.loadtxt(v18uixpath + 'LITbas_full.dat').astype(int)

basvs = {}
for bv in bas:
    try:
        basvs[bv[0]].append(bv[1])
    except:
        basvs[bv[0]] = [bv[1]]
sbas = sorted(basvs.items(), key=operator.itemgetter(0))

n3_inen_bdg(sbas, 1. / 2., costr, fn='INEN', pari=0, nzop=31, tni=11)

os.system(BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')
os.system('cp OUTPUT OUTPUT_full')

red_mod_3(
    typ='',
    max_coeff=11000,
    min_coeff=150,
    target_size=130,
    nbr_cycles=410,
    max_diff=0.001,
    ord=0,
    tniii=31,
    delpred=2,
    dr2executable=BINBDGpath + 'DR2END_AK.exe')

inen0 = [line for line in open('INEN')]

bas0 = []
nBV = int(inen0[7][4:8])
for bvn in range(1, 1 + nBV):
    bv = int(inen0[int(8 + 2 * (bvn - 1))][4:8])
    rel = np.argwhere(
        np.array(inen0[int(9 + 2 * (bvn - 1))].strip().split()).astype(
            int)).flatten() + 1
    bas0.append((bv, list(rel)))

with open(v18uixpath + 'LITbas_red.dat', 'w') as f:
    np.savetxt(f, [[jj[0], kk] for jj in bas0 for kk in jj[1]], fmt='%d')
f.close()

os.system(BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')
os.system('cp OUTPUT OUTPUT_red')

visbas(basispath=v18uixpath, widthpath=litpath3He, exepath=BINBDGpath)