from bridge import *
from three_particle_functions import *
from BasisVisualization import visbas
import operator

os.chdir(v18uixpath)
print('(working dir) %s' % v18uixpath)

th_low = 0.00001
th_up = .005
w_max = 1.5

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

intwLIT = [
    np.array(ln.split(';')).astype(float).tolist()
    for ln in open(litpath3He + 'intw3heLIT.dat')
]

relwLIT = [
    np.array(ln.split(';')).astype(float).tolist()
    for ln in open(litpath3He + 'relw3heLIT.dat')
]

iws = []
rws = []

bass = []

#print(' BV REL          wi          wr')
# loop through the basis
for bv in bas:
    for fr in range(len(intwLIT)):
        if bv[0] <= sum([len(ws) for ws in intwLIT[:fr + 1]]):
            iw = bv[0] - sum([len(ws) for ws in intwLIT[:fr]])
            iws.append(float(intwLIT[fr][iw - 1]))
            rws.append(float(relwLIT[fr][bv[1] - 1]))
            #print('%3d%3d%12.4f%12.4f' % (bv[0], bv[1], iws[-1], rws[-1]))
            # include bv in reduced basis if ``a criteria'' is satisfied
            if (((iws[-1] < th_up) | (rws[-1] < th_up)) &
                (w_max > iws[-1] > th_low) & (w_max > rws[-1] > th_low)):
                bass.append([bv[0], bv[1]])
            break

basvs = {}

for bv in bass:
    try:
        basvs[bv[0]].append(bv[1])
    except:
        basvs[bv[0]] = [bv[1]]

sbas = sorted(basvs.items(), key=operator.itemgetter(0))

with open(v18uixpath + 'LITbas_red.dat', 'w') as f:
    np.savetxt(f, [[jj[0], kk] for jj in sbas for kk in jj[1]], fmt='%d')
f.close()

n3_inen_bdg(sbas, 1. / 2., costr, fn='INEN', pari=0, nzop=31, tni=11)

os.system(BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')
os.system('cp OUTPUT OUTPUT_red')

visbas(basispath=v18uixpath, widthpath=litpath3He, exepath=BINBDGpath)