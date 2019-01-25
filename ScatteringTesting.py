#Scattering Testing

import numpy
import scipy
import IQHEDiag
import matplotlib
import matplotlib.pyplot as pyplot
import usefulTools

def findDispersion(spectrum, LMax):
    E = {0 : 0, 1 : 0}
    E_0 = spectrum[0][1]
    energyIndex = 2
    for L in range(2,LMax):
        E[L] = spectrum[energyIndex][1] - E_0
        energyIndex += len(usefulTools.generatePartitions(L))
    return E

def stateEnergyFromDispersion(p, E):
    return sum([E[L] for L in p])

def predictSpectrum(E, LMax, E_0):
    result = []
    for L in range(LMax):
        partitions = usefulTools.generatePartitions(L)
        result += [[L, E_0 + stateEnergyFromDispersion(p, E)] for p in partitions]
    return result

def spectrumCompareWithNoScatter(numericalSpectrum):
    LMax = numericalSpectrum[-1][0] + 1
    print(numericalSpectrum)
    E_0 = numericalSpectrum[0][1]
    E = findDispersion(numericalSpectrum, LMax)
    CFTSpectrum = predictSpectrum(E, LMax, E_0)
    print(CFTSpectrum)
    L1 = [item[0] for item in numericalSpectrum]
    E1 = [item[1] for item in numericalSpectrum]
    U = 0*(max(E1) - min(E1))/(LMax - 1)
    print("U = " + str(U))
    E1 = [item[1] + U*item[0] for item in numericalSpectrum]
    L2 = [item[0] for item in CFTSpectrum]
    E2 = [item[1] + U*item[0] for item in CFTSpectrum]
    Ls = [L for L in E if L > 0]
    EAbs = [-E[L] for L in E if L > 0]
    p = numpy.polyfit(Ls, EAbs, 4)
    for i in range(5):
        print(4-i)
        print(p[i])
    EPred = [numpy.polyval(p, L) for L in Ls]
    EDiff = [EAbs[i] - EPred[i] for i in range(len(Ls))]
    #pyplot.plot(Ls, EDiff)
    pyplot.plot(Ls, EPred, 'ko')
    pyplot.plot(Ls, EAbs)
    """
    pyplot.subplot(1,2,1)
    pyplot.xlabel("Angular momentum above ground state", fontsize = 24)
    pyplot.ylabel("Energy", fontsize = 24)
    pyplot.plot(L1, E1, 'bo')
    pyplot.hlines(E2, [i - 0.2 for i in L2], [i + 0.2 for i in L2])

    pyplot.subplot(1,2,2)
    pyplot.xlabel("L")
    pyplot.ylabel("|E(L)|")
    pyplot.plot(Ls, EAbs)
    """
    pyplot.show()

numericalSpectrum = IQHEDiag.findEnergiesForRangeOfL(40, 8, 1, 0)
spectrumCompareWithNoScatter(numericalSpectrum)
