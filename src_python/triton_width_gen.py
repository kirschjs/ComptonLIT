import numpy as np
import matplotlib.pyplot as plt
import os, re

import multiprocessing


def retrieve_he3_widths(inqua):

    relw = []
    intw = []
    frgm = []
    inq = [line for line in open(inqua)]

    lineNR = 0
    while lineNR < len(inq):
        if ((re.search('Z', inq[lineNR]) != None) |
            (re.search('z', inq[lineNR]) != None)):
            break
        lineNR += 1
    if lineNR == len(inq):
        print('no <Z> qualifier found in <INQUA>!')
        exit()

    while lineNR < len(inq):

        anziw = int(inq[lineNR].split()[0])
        anzbvLN = 2 * anziw if anziw <= 6 else 3 * anziw
        anzrw = int(inq[lineNR + 1].split()[1])
        anzrwLN = int(np.ceil(float(anzrw) / 6))

        frgm.append([anziw, anzrw])

        for iws in range(anziw):
            intw.append(float(inq[lineNR + 2 + iws].strip()))
        for rws in range(anzrwLN):
            for rw in inq[lineNR + 2 + int(inq[lineNR].split()[0]) +
                          rws].split():
                relw.append(rw)

        lineNR += anziw + anzrwLN + anzbvLN + 2

    iw = np.sort(np.array(intw))[::-1]
    rw = np.sort(np.array(relw).astype(float))[::-1]

    with open('intw3he.dat', 'wb') as f:
        for ws in iw:
            np.savetxt(f, [ws], fmt='%12.4f', delimiter=' ; ')
    f.close()
    with open('relw3he.dat', 'wb') as f:
        for ws in rw:
            np.savetxt(f, [ws], fmt='%12.4f', delimiter=' ; ')
    f.close()

    return iw, rw, frgm


def threebodywidths(zerl, means=[], nPERz=5):

    if means == []:
        # define centres for the means
        means = np.linspace(0.0, 2, nPERz)

    s = []
    for mu in means:
        s += [np.abs(np.random.normal(mu, 0.01, 400))]

    # for each cfg (e.g.: 101 - he3-no3)
    ws = []
    for zz in zerl:
        tmp = []
        for mean in s:
            tmp += [np.random.choice(mean, 1)]
        ws += [np.sort(np.array(tmp).flatten())[::-1]]

    return ws


def threebodywidths_int(zerl, means, nPERz=5):

    if means == []:
        # define centres for the means
        means = np.linspace(0.0, 2, nPERz)

    s = []
    for mu in means:
        s += [np.abs(np.random.normal(mu, 0.01, 400))]

    # for each cfg (e.g.: 101 - he3-no3)
    ws = []
    for zz in zerl:
        tmp = []
        for mean in means:
            tmp += [mean * (2 * np.random.rand() + 9) / 10]
        ws += [np.sort(np.random.choice(np.array(tmp).flatten(), nPERz))[::-1]]

    return ws


def perturbseedwidths(seeds, anz=20):

    if seeds == []:
        # define centres for the means
        means = np.linspace(0.0, 2, anz)

    # for each cfg (e.g.: 101 - he3-no3)
    ws = []
    while (len(ws) <= anz):
        for seed in seeds:
            ws.append(np.abs(np.random.normal(loc=seed, scale=0.4, size=None)))
            #ws.append(np.abs(seed * (2 * np.random.rand() + 9) / 10))

    return np.sort(ws)[::-1]


#count, bins, ignored = plt.hist(s, 100, normed=True, align='mid')
#x = np.linspace(min(bins), max(bins), 10000)
#pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2)) /
#       (x * sigma * np.sqrt(2 * np.pi)))
#plt.plot(x, pdf, linewidth=2, color='r')
#plt.axis('tight')
#plt.show()