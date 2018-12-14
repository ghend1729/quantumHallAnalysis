#CFT comparison
import math
import numpy
import scipy
from usefulTools import generatePartitions
import matplotlib
import matplotlib.pyplot as pyplot
import IQHEDiag

def zeroOrderEnergies(n, h_22, h_33):
    return 2*n*(n**2 - 1)*(h_33*(n**2 - 4) - h_22)

def stateEnergy(partition, h_22, h_33):
    return sum([zeroOrderEnergies(n, h_22, h_33) for n in partition])

def findParameters(spectrum):
    E_0 = spectrum[0][1]
    E_2 = spectrum[2][1]
    E_3 = spectrum[4][1]
    h_22 = (E_0 - E_2)/12
    h_33 = (E_3 - E_0 + 48*h_22)/240
    return E_0, h_22, h_33

def calcSpetrum(LMax, E_0, h_22, h_33):
    result = []
    for L in range(LMax):
        partitions = generatePartitions(L)
        result += [[L, E_0 + stateEnergy(p, h_22, h_33)] for p in partitions]
    return result

def spectrumCompare(numericalSpectrum):
    LMax = numericalSpectrum[-1][0] + 1
    print(numericalSpectrum)
    E_0, h_22, h_33 = findParameters(numericalSpectrum)
    CFTSpectrum = calcSpetrum(LMax, E_0, h_22, h_33)
    L1 = [item[0] for item in numericalSpectrum]
    E1 = [item[1] for item in numericalSpectrum]
    L2 = [item[0] for item in CFTSpectrum]
    E2 = [item[1] for item in CFTSpectrum]
    pyplot.xlabel("Delta L")
    pyplot.ylabel("E/(e^2/epsilon0/magnetic length/(4*pi))")
    pyplot.plot(L1, E1, 'bo')
    pyplot.plot(L2, E2, 'rx')
    pyplot.show()

spectrumCompare(IQHEDiag.findEnergiesForRangeOfL(30, 9, 1, 0))
