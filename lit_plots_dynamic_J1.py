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


rmin = 0
rmax = 10
suffix = 'sI15'

nl = 1
sigma = [s for s in open('sRange')]
sR = np.array(sigma[1].split(';')).astype(float)
ns = len(sR)
sI = float(sigma[0])
lit_wrels = [float(li) for li in open('wLIT.dat')]
nw = len(lit_wrels)

litcoeffs10 = np.array([
    complex(float(litc.split(';')[0]), float(litc.split(';')[1]))
    for litc in open("COEFFS_LIT_J1mJ0")
])

print('|Re[sigma]| = %d' % ns)
print('|w(LIT)|    = %d' % nw)
print('|COEFFs|    = %d' % len(litcoeffs10))
litcoeffs10 = np.reshape(litcoeffs10, (len(sR), -1))

rspace = np.linspace(rmin, rmax, 100)

# choose energies based on the behavior of the determinant (see cpts_..._RGM.nb)

wfktJ1 = [[psiL(rr, 1, lit_wrels, np.real(litcoeff)) for rr in rspace]
          for litcoeff in litcoeffs10]

f = plt.figure(figsize=(10, 6))
f.suptitle(
    r'$\Psi_{LIT}^{1^-}\;\;\;\;\sigma_i=%2.1f$ MeV $\;\;\;\;\sigma_r\in[%2.1f,%2.1f]$ MeV  (dark$\to$light)'
    % (sI, float(sR[0]), float(sR[-1])),
    fontsize=14)
ax1 = f.add_subplot(111)

ax1.set_xlabel(r'$r$ [fm]', fontsize=12)
ax1.set_ylabel(r'$\psi_{1^-}(r)$ [MeV$^{-3/2}$]', fontsize=12)
colormap = cm.Greys(np.linspace(0, 1, len(wfktJ1)))[::-1]

[
    ax1.plot(
        rspace,
        wfktJ1[wfkt],
        label='$\sigma_r=%4.2f, n=%d$' % (float(sR[wfkt]), nw),
        c=colormap[wfkt]) for wfkt in range(len(wfktJ1))
]

line, = ax1.plot(rspace, wfktJ1[0], 'r-')

for i in range(len(wfktJ1))[:-1]:
    #xdata.append(rspace[i])
    #ydata.append(ysample[i])
    #print(xdata, ydata)
    #line.set_title(r'$\sigma_r=%4.2f MeV$' % i)
    line.set_xdata(rspace)
    line.set_ydata(wfktJ1[i])
    if int(i) % 3 == 0:
        ax1.annotate(
            r'$\sigma_r = %d$' % i,
            xy=(0.94, 1 - i / len(wfktJ1)),
            xycoords='figure fraction')
    plt.draw()
    plt.pause(1e-7)
    time.sleep(.9)

ax1.yaxis.set_major_formatter(FormatStrFormatter('%1.0e'))
#plt.legend(loc='center left', numpoints=1, fontsize=18,
#           bbox_to_anchor=(1.1, .4))
strFile = 'LIT_waveFkt_J1_%s.pdf' % suffix
if os.path.isfile(strFile):
    os.remove(strFile)
plt.savefig(strFile)
plt.savefig(strFile)

plt.show()