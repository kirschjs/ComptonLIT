import numpy as np
import matplotlib.pyplot as plt
import os

import multiprocessing


def threebodywidths(zerl, nPERz=5):
    # define centres for the means
    means = np.linspace(0.1, 2, nPERz)
    s = []
    for mu in means:
        s += [np.random.lognormal(mu, 0.5, 400)]

    # for each cfg (e.g.: 101 - he3-no3)
    ws = []
    for zz in zerl:
        tmp = []
        for mean in s:
            tmp += [np.random.choice(mean, 1)]
        ws += [np.sort(np.array(tmp).flatten())[::-1]]

    return ws


#count, bins, ignored = plt.hist(s, 100, normed=True, align='mid')
#x = np.linspace(min(bins), max(bins), 10000)
#pdf = (np.exp(-(np.log(x) - mu)**2 / (2 * sigma**2)) /
#       (x * sigma * np.sqrt(2 * np.pi)))
#plt.plot(x, pdf, linewidth=2, color='r')
#plt.axis('tight')
#plt.show()