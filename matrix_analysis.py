import os, sys

import numpy as np
from scipy import linalg

infil = '/home_th/kirscher/kette_repo/source/ComptonLIT/av18_deuteron/norm-ham-litME-0-'

normham = np.array([float(line) for line in open(infil)])

EN = np.reshape(normham[1:int(normham[0]**2 + 1)], (int(normham[0]), -1))
H = np.reshape(normham[int(normham[0]**2 + 1):], (int(normham[0]), -1))

evg, evrg = linalg.eig(H, EN)
ev, evr = linalg.eig(H)
ev2, evr2 = linalg.eig(EN)

print('NORM EV:\n', np.sort(ev2)[::-1])
print('HAMI EV:\n', np.sort(ev)[::-1])
print('(HN) EV:\n', np.sort(evg)[::-1])