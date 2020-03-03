from scipy.special import binom
import numpy as np


def clebgor(j1, j2, j, m1, m2, m):
    """
    Parameters:    j1, j2, j: the angular momenta input
                   m1, m2, m: the z components of angular momenta input
    Returns:       The numerical value of the Clebsch-Gordan coeffcient.
    Remarks:       Note that in the sum none of the binomial coeffcients
                   can have negative values.  Thus, zmin is there to make
                   sure that the sums have a cut-off.
    """
    if ((abs(j1 - j2) > j) | (abs(j1 + j2) < j) | (m1 + m2 != m)):
        return
    zmin = int(min([j1 - m1, j2 + m2]))
    J = j1 + j2 + j
    return (int(m1 + m2 == m) * int(np.abs(j1 - j2) <= j <= (j1 + j2)) *
            int(np.abs(m1) <= j1) * int(np.abs(m2) <= j2) *
            int(np.abs(m) <= j) * int((j1 + m1) >= 0.0) * int((j2 + m2) >= 0.0)
            * int((j + m) >= 0.0) * int(J >= 0) * np.sqrt(
                binom(2 * j1, J - 2 * j) * binom(2 * j2, J - 2 * j) /
                (binom(J + 1, J - 2 * j) * binom(2 * j1, j1 - m1) * binom(
                    2 * j2, j2 - m2) * binom(2 * j, j - m))) *
            np.sum([(-1)**z * binom(J - 2 * j, z) * binom(
                J - 2 * j2, j1 - m1 - z) * binom(J - 2 * j1, j2 + m2 - z)
                    for z in range(zmin + 1)]))