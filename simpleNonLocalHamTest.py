#Simple non-local hamiltonain Testing
"""
We try including and B-O term to try and account for the linear dispersion
seen in the coulomb interaction.
"""

import numpy
import scipy
import matplotlib
import matplotlib.pyplot as pyplot
import IQHEDiag
import FQHEDiag
from usefulTools import generatePartitions
from numpy.linalg import inv

def extractParameters(spectrum):
    E_0 = spectrum[0][1]
    E_2 = spectrum[2][1]
    E_3 = spectrum[4][1]
    E_4 = spectrum[7][1]
    energyArray = [E_2 - E_0, E_3 - E_0, E_4 - E_0]
    r22 = [-2*n*(n**2 - 1) for n in range(2, 5)]
    r33 = [2*n*(n**2 - 4)*(n**2 - 1) for n in range(2, 5)]
    rg = [n*(n - 1) for n in range(2, 5)]
    fitterMatrix = numpy.transpose([r22, r33, rg])
    coeficients = inv(fitterMatrix).dot(energyArray)
    h_22 = coeficients[0]
    h_33 = coeficients[1]
    g = coeficients[2]
    return E_0, g, h_22, h_33

def energySingleParticle(n, g, h_22, h_33):
    return g*n*(n-1) - h_22*2*n*(n**2 - 1) + h_33*2*n*(n**2 - 1)*(n**2 - 4)

def stateEnergy(p, g, h_22, h_33):
    return sum([energySingleParticle(n, g, h_22, h_33) for n in p])

def calcPredictedSpectrum(LMax, g, h_22, h_33, E_0):
    predSpectrum = []
    for L in range(LMax):
        partitions = generatePartitions(L)
        predSpectrum += [[L, E_0 + stateEnergy(p, g, h_22, h_33)] for p in partitions]
    return predSpectrum

def spectrumCompare(numericalSpectrum):
    LMax = numericalSpectrum[-1][0] + 1
    print(numericalSpectrum)
    E_0, g, h_22, h_33 = extractParameters(numericalSpectrum)
    CFTSpectrum = calcPredictedSpectrum(LMax, g, h_22, h_33, E_0)
    print(CFTSpectrum)
    L1 = [item[0] for item in numericalSpectrum]
    E1 = [item[1] for item in numericalSpectrum]
    U = 0*(max(E1) - min(E1))/(LMax - 1)
    print("U = " + str(U))
    E1 = [item[1] + U*item[0] for item in numericalSpectrum]
    L2 = [item[0] for item in CFTSpectrum]
    E2 = [item[1] + U*item[0] for item in CFTSpectrum]
    pyplot.xlabel("Angular momentum above ground state", fontsize = 24)
    pyplot.ylabel("Energy/w_0", fontsize = 24)
    pyplot.plot(L1, E1, 'bo')
    pyplot.hlines(E2, [i - 0.2 for i in L2], [i + 0.2 for i in L2])
    pyplot.show()

numericalSpectrum = IQHEDiag.findEnergiesForRangeOfL(50, 8  , 1, 0)
spectrumCompare(numericalSpectrum)
