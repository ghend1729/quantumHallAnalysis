#Scaling Testing
"""
This script is to test how well non-local terms have been acounted for by
different proposed hamiltonians by looking at the scaling.
"""

import numpy
import math
import IQHEDiag
import scipy
from CfTComparison import *
import matplotlib
import matplotlib.pyplot as pyplot

def generateHs():
    h_22NormedVals = []
    h_33NormedVals = []
    m = 1
    for N in range(40, 91):
        spectrum = IQHEDiag.findEnergiesForRangeOfL(N, 5, 1, 0)
        E_0 = spectrum[0][1]
        E_2 = spectrum[2][1]
        E_3 = spectrum[4][1]
        E_4 = spectrum[7][1]
        energyArray = [E_2 - E_0, E_3 - E_0, E_4 - E_0]
        r22 = [-2*n*(n**2 - 1) for n in range(2, 5)]
        r33 = [2*n*(n**2 - 4)*(n**2 - 1) for n in range(2, 5)]
        rXi = [xi(n, N, m) for n in range(2, 5)]
        fitterMatrix = numpy.transpose([r22, r33, rXi])
        coeficients = inv(fitterMatrix).dot(energyArray)
        h_22 = coeficients[0]
        h_33 = coeficients[1]
        f = coeficients[2]
        h_22NormedVals.append(h_22*math.sqrt(N)**3)
        h_33NormedVals.append(h_33*math.sqrt(N)**5)
    return h_22NormedVals, h_33NormedVals

hval2, hval3 = generateHs()

x = [i for i in range(40, 91)]

pyplot.subplot(2,1,1)
pyplot.xlabel("N", fontsize = 18)
pyplot.ylabel("h_22/sqrt(N)^(-3)", fontsize = 18)
pyplot.plot(x, hval2, 'ko')

pyplot.subplot(2,1,2)
pyplot.xlabel("N", fontsize = 18)
pyplot.ylabel("h_33/sqrt(N)^(-5)", fontsize = 18)
pyplot.plot(x, hval3, 'ko')

pyplot.show()
