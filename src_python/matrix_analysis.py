import os, sys

import numpy as np
from scipy import linalg

home = os.getenv("HOME")

infil = home + '/kette_repo/ComptonLIT/av18_deuteron/norm-ham-litME-0-'

normham = np.array([float(line) for line in open(infil)])

Norm = np.transpose(
    np.reshape(normham[1:int(normham[0]**2 + 1)], (int(normham[0]), -1)))
Ham = np.reshape(normham[int(normham[0]**2 + 1):], (int(normham[0]), -1))

ev_gen, evr_gen = linalg.eig(Ham, Norm)
ev_h, evr_h = linalg.eig(Ham)
ev_n, evr_n = linalg.eig(Norm)

print(Norm)
print(Ham)

#print('NORM     EV:\n', np.sort(ev_n))
#print('H        EV:\n', np.sort(ev_h))
print('(Hv=lNv) EV:\n', np.sort(ev_gen))
exit()
print('\n', evr_gen)