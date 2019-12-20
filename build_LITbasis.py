from bridge import *
from three_particle_functions import *
from triton_width_gen import *
import operator
from BasisVisualization import visbas

os.chdir(v18uixpath)

min_eucl_pair_dist = 1.5

intwLIT = [
    np.array(ln.split(';')).astype(float).tolist()
    for ln in open(litpath3He + 'intw3heLIT.dat')
]

relwLIT = [
    np.array(ln.split(';')).astype(float).tolist()
    for ln in open(litpath3He + 'relw3heLIT.dat')
]

bas = []
for nz in range(len(intwLIT)):
    allindices = [[w1, w2] for w1 in range(len(intwLIT[nz]))
                  for w2 in range(len(relwLIT[nz]))]
    e0 = np.random.choice(range(len(allindices)))
    set0 = [allindices[e0]]
    allindices.remove(allindices[e0])
    while len(allindices) > 0:
        e0 = np.random.choice(range(len(allindices)))
        reject = False
        for ei in range(len(set0)):
            w1 = [
                intwLIT[nz][allindices[e0][0]], relwLIT[nz][allindices[e0][1]]
            ]
            w2 = [intwLIT[nz][set0[ei][0]], relwLIT[nz][set0[ei][1]]]
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

    offset_int = sum([len(iws) for iws in intwLIT[:nz]])
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

visbas(basispath=v18uixpath, widthpath=litpath3He, exepath=BINBDGpath)

n3_inen_bdg(sbas, 1. / 2., costr, fn='INEN', pari=0, nzop=31, tni=11)
os.system(BINBDGpath + 'DR2END_AK.exe && grep -A 3 \'EIGENWER\' OUTPUT')
