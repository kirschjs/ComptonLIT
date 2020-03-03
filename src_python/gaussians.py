import os
import numpy as np
import matplotlib.pyplot as plt
from parameters_and_constants import *


def sparsify(menge, mindist):

    lockereMenge = []
    menge = np.sort(menge)
    print(menge, '\n')
    npointer = 1
    nref = 0

    while (npointer < len(menge)):

        if (np.abs(float(menge[nref]) - float(menge[npointer])) > mindist):
            lockereMenge.append(float(menge[nref]))
            nref = npointer
            npointer = nref + 1

        else:
            npointer += 1

    if (np.abs(float(menge[-1]) - float(lockereMenge[-1])) > mindist):
        lockereMenge.append(float(menge[-1]))

    return np.sort(lockereMenge)


anz = 15

loc, scale = 1., 2.
wLAPLACE = np.sort(np.abs(np.random.laplace(loc, scale, anz)))

print(sparsify(wLAPLACE, 0.1))
exit()

n_bins = 30
maxvalue = 13.0

wLEGE = np.unique(
    np.sort(np.abs(np.real(np.polynomial.legendre.leggauss(anz)[0]))))
wLEGEscaled = wLEGE * maxvalue / np.max(wLEGE)
wHERM = np.unique(
    np.sort(np.abs(np.real(np.polynomial.hermite.hermgauss(anz)[0]))))
wHERMscaled = wHERM * maxvalue / np.max(wHERM)
wLAGU = np.unique(
    np.sort(np.abs(np.real(np.polynomial.laguerre.laggauss(anz)[0]))))
wLAGUscaled = wLAGU * maxvalue / np.max(wLAGU)
wCHEB = np.unique(
    np.sort(np.abs(np.real(np.polynomial.chebyshev.chebgauss(anz)[0]))))
wCHEBscaled = wCHEB * maxvalue / np.max(wCHEB)

print(wCHEB, wHERM)

fig, ax = plt.subplots(1, 2, sharey=True, tight_layout=True)

ax[0].hist(wLEGEscaled, bins=n_bins, label=r'Legendre')
ax[0].hist(wLAGUscaled, bins=n_bins, label=r'Lagurre')
ax[0].hist(wCHEBscaled, bins=n_bins, label=r'Chebyshev')
ax[1].hist(w12, bins=n_bins, label=r'w12')

ax[0].legend(loc='best')
plt.show()