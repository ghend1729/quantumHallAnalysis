#Dispersion/Scattering Testing

import numpy
import scipy
import math
import IQHEDiag
import FQHEDiag
import matplotlib
import matplotlib.pyplot as pyplot
import usefulTools
import scipy.optimize
import pickle

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

def delta(x):
    if x == 1:
        return 1
    else:
        return 0

def f(n, g, h, a, b):
    return g*(n-1)*n + h*(n-1) + a*n*(n-1)*(n-2)

def spectrumCompareWithNoScatter(numericalSpectrum, N):
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
    LFit = [L for L in E if L > 0]
    EFit = [-E[L] for L in E if L > 0]
    g, h, a, b = scipy.optimize.curve_fit(f, LFit, EFit)[0]
    print(g, h, a, b)
    EPred = [f(L, g, h, a, b) for L in Ls]
    EDiff = [(EAbs[i] - EPred[i])/max(EAbs) for i in range(len(EAbs))]
    EDiff2 = [EAbs[i] - EAbs[i-1] for i in range(1, len(Ls))]
    LDiff2 = [Ls[i] for i in range(1, len(Ls))]
    findPeak(LDiff2, EDiff2)
    pyplot.title("N = 339 Difference Graph", fontsize = 18)
    pyplot.xlabel("N", fontsize = 20)
    pyplot.ylabel("|E(n) - E(n-1)|", fontsize = 20)
    pyplot.plot(LDiff2, EDiff2, 'ko')
    #pyplot.plot(Ls, EPred)
    """
    pyplot.subplot(1,2,1)
    pyplot.xlabel("Angular momentum above ground state", fontsize = 24)
    pyplot.ylabel("Energy", fontsize = 24)
    pyplot.plot(L1, E1, 'bo')
    pyplot.hlines(E2, [i - 0.2 for i in L2], [i + 0.2 for i in L2])

    pyplot.subplot(1,2,2)
    pyplot.xlabel("L")
    pyplot.ylabel("|E(L)|")
    pyplot.plot(Ls, EAbs, 'ko')
    pyplot.plot(Ls, EPred)
    """
    pyplot.show()

def pow(x, a, b):
    return b*x**a

def scaleTest(spectra):
    Ns = []
    hsScaled = []
    gsScaled = []
    hs = []
    gs = []
    As = []
    asScaled = []
    for N in range(150, 200):
        spectrum = spectra[(N, 11)]
        E = findDispersion(spectrum, 11)
        Ls = [L for L in E if L > 1 and L < 7]
        EAbs = [-E[L] for L in E if L > 1 and L < 7]
        g, h, a, b = scipy.optimize.curve_fit(f, Ls, EAbs)[0]
        hsScaled.append(h*math.sqrt(N))
        gsScaled.append(g*N**(1/2))
        Ns.append(N)
        hs.append(h)
        gs.append(g)
        As.append(a)
        asScaled.append(a*N**(2/3))

    a, b = scipy.optimize.curve_fit(pow, Ns, As)[0]
    print(a, b)

    pyplot.subplot(1,2,1)
    pyplot.xlabel("N", fontsize = 20)
    pyplot.ylabel("g")
    pyplot.plot(Ns, gsScaled)

    pyplot.subplot(1,2,2)
    pyplot.xlabel("N", fontsize = 20)
    pyplot.ylabel("h")
    pyplot.plot(Ns, hsScaled)

    pyplot.show()

def generateData():
    spectra = {}
    #first set for scaling testing
    for N in range(70, 200):
        spectrum = IQHEDiag.findEnergiesForRangeOfL(N, 11, 1, 0)
        spectra[(N, 11)] = spectrum
    #Second set for fitting testing
    Ns = [30, 70, 140]
    for N in Ns:
        spectrum = IQHEDiag.findEnergiesForRangeOfL(N, 19, 1, 0)
        spectra[(N, 19)] = spectrum
    #Now save the spectra
    spectraFile = open("BigSpectraCollection.p", 'wb')
    pickle.dump(spectra, spectraFile)
    spectraFile.close()

def findPeak(LDiff, EDiff):
    approxPeakIndex = EDiff.index(max(EDiff))
    fitDataL = LDiff[approxPeakIndex - 1: approxPeakIndex + 2]
    fitDataEDiff = EDiff[approxPeakIndex - 1: approxPeakIndex + 2]
    p = numpy.polyfit(fitDataL, fitDataEDiff, 2)
    return -p[1]/(2*p[0])

def peakAnalysis():
    Ns = []
    LPeaks = []
    for N in range(30, 140):
        spectrum = spectra[(N, 15)]
        E = findDispersion(spectrum, 15)
        Ls = [L for L in E if L > 1]
        EAbs = [-E[L] for L in E if L > 1]
        EDiff2 = [EAbs[i] - EAbs[i-1] for i in range(1, len(Ls))]
        LDiff2 = [Ls[i] for i in range(1, len(Ls))]
        peak = findPeak(LDiff2, EDiff2)
        LPeaks.append(peak)
        Ns.append(N)
    a, b = scipy.optimize.curve_fit(pow, Ns, LPeaks)[0]
    print(a, b)
    pyplot.ylabel("L value of peak", fontsize = 18)
    pyplot.xlabel("N", fontsize = 18)
    pyplot.plot(Ns, LPeaks)
    pyplot.show()

spectraFile = open("BigSpectraCollection.p", 'rb')
spectra = pickle.load(spectraFile)
spectraFile.close()
spectrumCompareWithNoScatter(spectra[(339, 11)], 339)
#scaleTest(spectra)
#peakAnalysis()
"""
for N in range(200, 301):
    spectrum = IQHEDiag.findEnergiesForRangeOfL(N, 11, 1, 0)
    spectra[(N, 11)] = spectrum
IQHEDiag.dumpRequest()
spectraFile2 = open("BigSpectraCollection.p", 'wb')
pickle.dump(spectra, spectraFile2)
spectraFile2.close()
"""
