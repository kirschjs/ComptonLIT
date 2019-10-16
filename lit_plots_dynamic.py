import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colorbar, colors
from matplotlib.ticker import FormatStrFormatter
from matplotlib.collections import LineCollection
from bridge import *
import time
import random

os.chdir(av18path)


def psiL(r, l, widths, coeffs):
    psi = np.array([
        coeffs[n] * r**l * np.exp(-widths[n] * r**2)
        for n in range(min(len(coeffs), len(widths)))
    ])
    return np.sum(psi)


deut_expansion = [float(li) for li in open('COEFF')]
clean_wrels = [float(li) for li in open('w0.dat')]
print('dim(deuteron) = ', len(deut_expansion))
print('anz(w0_clean) = ', len(clean_wrels))

rmin = 0
rmax = 10
suffix = 'sI15'

inen = [line for line in open('inen_b')]
Sindices = np.nonzero(np.array(inen[6].split()).astype(int))
Dindices = np.nonzero(np.array(inen[8].split()).astype(int))
Swidths = np.take(clean_wrels, Sindices)[0]
Dwidths = np.take(clean_wrels, Dindices)[0]
print('S-wave dim(deuteron) = ', len(Sindices[0]))
print('D-wave dim(deuteron) = ', len(Dindices[0]))

nl = 1
sigma = [s for s in open('sRange')]
sR = np.array(sigma[1].split(';')).astype(float)
ns = len(sR)
sI = float(sigma[0])
lit_wrels = [float(li) for li in open('wLIT.dat')]
nw = len(lit_wrels)
litcoeffs00 = np.array([
    complex(float(litc.split(';')[0]), float(litc.split(';')[1]))
    for litc in open("COEFFS_LIT_J0mJ0")
])
print('|Re[sigma]| = %d' % ns)
print('|w(LIT)|    = %d' % nw)
print('|COEFFs|    = %d' % len(litcoeffs00))
litcoeffs00 = np.reshape(litcoeffs00, (len(sR), -1))

rspace = np.linspace(rmin, rmax, 100)

# choose energies based on the behavior of the determinant (see cpts_..._RGM.nb)

wfktJ0 = [[psiL(rr, 1, lit_wrels, np.real(litcoeff)) for rr in rspace]
          for litcoeff in litcoeffs00]

f = plt.figure(figsize=(10, 6))
f.suptitle(
    r'$\Psi_{LIT}^{0^-}\;\;\;\;\sigma_i=%2.1f$ MeV $\;\;\;\;\sigma_r\in[%2.1f,%2.1f]$ MeV (dark$\to$light)'
    % (sI, float(sR[0]), float(sR[-1])),
    fontsize=14)

ax1 = f.add_subplot(111)
colormap = cm.Greys(np.linspace(0, 1, len(wfktJ0)))[::-1]

ax1.set_xlabel(r'$r$ [fm]', fontsize=12)
ax1.set_ylabel(r'$\psi_{0^-}(r)$ [MeV$^{-3/2}$]', fontsize=12)

[
    ax1.plot(
        rspace,
        wfktJ0[wfkt],
        label='$\sigma_r=%4.2f, n=%d$' % (float(sR[wfkt]), nw),
        c=colormap[wfkt]) for wfkt in range(len(wfktJ0))
]

line, = ax1.plot(rspace, wfktJ0[0], 'r-')

energysel = [1, 4, 23, 45, 76, 90]

for i in range(len(wfktJ0)):
    #xdata.append(rspace[i])
    #ydata.append(ysample[i])
    #print(xdata, ydata)
    #line.set_title(r'$\sigma_r=%4.2f MeV$' % i)
    line.set_xdata(rspace)
    line.set_ydata(wfktJ0[i])
    if int(i) % 3 == 0:
        ax1.annotate(
            r'$\sigma_r = %d$' % i,
            xy=(0.94, 1 - i / len(wfktJ0)),
            xycoords='figure fraction')
    plt.draw()
    plt.pause(1e-7)
    time.sleep(0.1)

ax1.yaxis.set_major_formatter(FormatStrFormatter('%1.0e'))
#plt.legend(loc='center left', numpoints=1, fontsize=18,
#           bbox_to_anchor=(1.1, .4))
strFile = 'LIT_waveFkt_J0_%s.pdf' % suffix
if os.path.isfile(strFile):
    os.remove(strFile)
plt.savefig(strFile)
plt.savefig(strFile)

plt.show()