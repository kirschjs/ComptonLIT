import os, sys

import numpy as np
from scipy import linalg

infil = '/home_th/kirscher/kette_repo/source/ComptonLIT/av18_deuteron/norm-ham-litME-0-'

normham = np.array([float(line) for line in open(infil)])

EN = np.reshape(normham[1:int(normham[0]**2 + 1)], (int(normham[0]), -1))
H = np.reshape(normham[int(normham[0]**2 + 1):], (int(normham[0]), -1))

ev, evr = linalg.eig(H, EN)

print(np.sort(ev))