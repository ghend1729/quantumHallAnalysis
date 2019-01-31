#Scattering Testing

import numpy
import scipy
import math
import IQHEDiag
import FQHEDiag
import matplotlib
import matplotlib.pyplot as pyplot
import usefulTools
import scipy.optimize

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

def f(n, g, h, h_22):
    return g*n*(n-1) + h*(n-1) + h_22*n*(n**2 - 1)

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
    g, h, h_22 = scipy.optimize.curve_fit(f, Ls, EAbs)[0]
    print(g, h, h_22)
    EPred = [f(L, g, h, h_22) for L in Ls]
    EDiff = [EAbs[i] - EPred[i] for i in range(len(Ls))]
    pyplot.xlabel("L", fontsize = 20)
    pyplot.ylabel("|E(L)|", fontsize = 20)
    pyplot.plot(Ls, EAbs, 'ko')
    pyplot.plot(Ls, EPred)
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

def scaleTest():
    Ns = []
    hsScaled = []
    h_22sScaled = []

    for N in range(70, 140):
        spectrum = IQHEDiag.findEnergiesForRangeOfL(N, 11, 1, 0)
        E = findDispersion(spectrum, 11)
        Ls = [L for L in E if L > 0]
        EAbs = [-E[L] for L in E if L > 0]
        g, h, h_22 = scipy.optimize.curve_fit(f, Ls, EAbs)[0]
        hsScaled.append(h*(math.sqrt(N)))
        h_22sScaled.append(h_22*(math.sqrt(N)**3))
        Ns.append(N)

    pyplot.subplot(1,2,1)
    pyplot.xlabel("N", fontsize = 20)
    pyplot.ylabel("h/sqrt(N)^-1")
    pyplot.plot(Ns, hsScaled)

    pyplot.subplot(1,2,2)
    pyplot.xlabel("N", fontsize = 20)
    pyplot.ylabel("h_22/sqrt(N)^-3")
    pyplot.plot(Ns, h_22sScaled)

    pyplot.show()

#spectrum = IQHEDiag.findEnergiesForRangeOfL(200, 11, 100, 0)
#spectrumCompareWithNoScatter(spectrum)
scaleTest()
