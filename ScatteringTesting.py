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
    print(E)
    result = []
    for L in range(LMax):
        partitions = usefulTools.generatePartitions(L)
        result += [[L, E_0 + stateEnergyFromDispersion(p, E)] for p in partitions]
    return result

def averageError(actualSpec, predSpec):
    count = 0
    total = 0
    E_0 = actualSpec[0][1]
    deltaE = E_0 - min([x[1] for x in actualSpec])
    LMax = max([x[0] for x in actualSpec]) + 1
    for L in range(LMax):
        energies = sorted([x[1] for x in actualSpec if x[0] == L])[:-1]
        for E in energies:
            total += min([abs(E - x[1]) for x in predSpec if x[0] == L])/deltaE
            count += 1
    return 100*total/count

def delta(x):
    if x == 1:
        return 1
    else:
        return 0

def f(n, g, h, a, b):
    return (h + g*n + a*(n**3) + b*(n**5))*(n-1)

def spectrumCompareWithNoScatter(numericalSpectrum, N):
    LMax = numericalSpectrum[-1][0] + 1
    print(numericalSpectrum)
    E_0 = numericalSpectrum[0][1]
    E = findDispersion(numericalSpectrum, LMax)
    CFTSpectrum = predictSpectrum(E, LMax, E_0)
    print(CFTSpectrum)
    L1 = [item[0] for item in numericalSpectrum if item[0] < 8]
    E1 = [item[1] for item in numericalSpectrum if item[0] < 8]
    U = 0*(max(E1) - min(E1))/(LMax - 1)
    print("U = " + str(U))
    E1 = [item[1] + U*item[0] for item in numericalSpectrum if item[0] < 8]
    L2 = [item[0] for item in CFTSpectrum if item[0] < 8]
    E2 = [item[1] + U*item[0] for item in CFTSpectrum if item[0] < 8]
    Ls = [L for L in E if L > 0]
    EAbs = [-E[L] for L in E if L > 0]
    LFit = [L for L in E if L > 0]
    EFit = [-E[L] for L in E if L > 0]
    
    EDiff2 = [EAbs[i] - EAbs[i-1] for i in range(1, len(Ls))]
    LDiff2 = [Ls[i] for i in range(1, len(Ls))]

    print("")
    for x in numericalSpectrum:
        print([x[0], E_0 - x[1]])
    """
    Ns = []
    errorsList = []

    for N2 in range(30, 340, 20):
        spectrum = spectra[(N2, 11)]
        LMax = spectrum[-1][0] + 1
        E_0 = spectrum[0][1]
        E = findDispersion(spectrum, LMax)
        CFTSpectrum = predictSpectrum(E, LMax, E_0)
        Ns.append(N2)
        errorsList.append(averageError(spectrum, CFTSpectrum))

    (a, b), errorMatrix = scipy.optimize.curve_fit(pow, Ns, errorsList)
    aError = math.sqrt(errorMatrix[0][0])
    bError = math.sqrt(errorMatrix[1][1])
    errorListPred = [pow(x, a, b) for x in range(30, 339)]
    NsPred = [x for x in range(30, 339)]
    print(a, b)
    print(aError, bError)
    """
    ax = pyplot.subplot(121)
    ax.tick_params(labelsize = 15)
    pyplot.title("$N = " + str(N) + "$\alpha = 0.1$ Free Boson Test", fontsize = 22)
    pyplot.xlabel("$\Delta L$", fontsize = 20)
    pyplot.ylabel("$\Delta E/(q^2/4\pi\epsilon_0l_B)$", fontsize = 20)
    pyplot.plot(L1, E1, 'ko', label = "$\delta\hat{H}$")
    pyplot.hlines(E2, [i - 0.3 for i in L2], [i + 0.3 for i in L2], label = "$\hat{H}$")
    pyplot.text(0.98, 0.02, "$U_0 = " + str(round(U/2, 4)) + "q^2/4\pi\epsilon_0l_B^3$", horizontalalignment='right', fontsize = 20, transform=ax.transAxes)
    pyplot.legend(fontsize = 22)
    
    ax2 = pyplot.subplot(122)
    ax2.tick_params(labelsize=15)
    pyplot.title("Mode energy vs n", fontsize=25)
    pyplot.xlabel("n", fontsize=22)
    pyplot.ylabel("E", fontsize = 22)
    pyplot.plot(LDiff2, EDiff2, label = "ACTUAL ERROR")
    
    pyplot.show()

def pow(x, a, b):
    return b*x**a

def bn(n, p0, p1, p2):
    return (n-1)*(p0 + p1*n + p2*n**2)

def scaleTest(spectra, f):
    Ns = []
    fs = [[] for i in range(9)]
    for N in range(150, 340):
        spectrum = spectra[(N, 11)]
        E = findDispersion(spectrum, 11)
        Ls = [L for L in E if L > 0]
        EAbs = [-E[L] for L in E if L > 0]
        Ns.append(N)
        for i in range(9):
            fs[i].append(EAbs[i+1])

    bs = []
    a2 = []
    bErrors = []
    aErrors = []
    for i in range(9):
        (a, b), M = scipy.optimize.curve_fit(pow, Ns, fs[i])
        print(a,b)
        bs.append(b)
        a2.append(a)
        bErrors.append(math.sqrt(M[1][1]))
        aErrors.append(math.sqrt(M[0][0]))

    ns = [i+2 for i in range(9)]
    for i in range(9):
        print(bErrors[i], aErrors[i])

    ax = pyplot.subplot(1,2,1)
    ax.tick_params(labelsize = 15)
    pyplot.title("m=1 g(2) vs N", fontsize = 22)
    pyplot.xlabel("N", fontsize = 20)
    pyplot.ylabel("$g(2,N)$", fontsize = 20)
    pyplot.plot([Ns[i] for i in range(0, len(Ns), 15)], [fs[0][i] for i in range(0, len(fs[0]), 15)], 'ko', label = "NUMERICAL")
    pyplot.plot(Ns, [pow(x, a2[0], bs[0]) for x in Ns], label = "$g(2,N) = b(2)N^{a(2)}$ FIT")
    pyplot.text(0.98, 0.75, "a(2) = " + str(round(a2[0], 4)) + "\n" + "b(2) = " + str(round(bs[0], 3)), horizontalalignment='right', verticalalignment='center', fontsize=20, transform=ax.transAxes)
    pyplot.legend(fontsize = 20)

    ns = [1] + ns
    bs = [0] + bs

    (p0, p1, p2), M = scipy.optimize.curve_fit(bn, ns, bs)
    print("")
    print(p0, math.sqrt(M[0][0]))
    print(p1, math.sqrt(M[1][1]))
    print(p2, math.sqrt(M[2][2]))
    print("")
    ax2 = pyplot.subplot(1,2,2)
    ax2.tick_params(labelsize = 15)
    pyplot.title("m = 1 b(n) Data & Fit", fontsize = 22)
    pyplot.xlabel("n", fontsize = 20)
    pyplot.ylabel("$b(n)$", fontsize = 20)
    pyplot.plot(ns, bs, 'ko', label = "NUMERICAL")
    pyplot.plot(ns, [bn(n, p0, p1, p2) for n in ns], label = "POLYNOMIAL FIT")
    pyplot.legend(fontsize = 20)

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

def differenceGraphPlotter():
    differenceData = []
    Ns = [30, 45, 60, 90, 120]
    for N in Ns:
        spectrum = spectra[(N, 15)]
        E = findDispersion(spectrum, 15)
        Ls = [L for L in E if L > 0]
        EAbs = [-E[L] for L in E if L > 0]
        EDiff2 = [math.sqrt(N)*(EAbs[i] - EAbs[i-1]) for i in range(1, len(Ls))]
        LDiff2 = [Ls[i] for i in range(1, len(Ls))]
        differenceData.append(EDiff2)

    pyplot.xlim((2, 16.7))
    pyplot.title("$h(n)$ for a range of $N$", fontsize = 20)
    pyplot.tick_params(labelsize = 12)
    pyplot.xlabel("$n$", fontsize = 18)
    pyplot.ylabel("$h(n)$", fontsize = 18)
    for i in range(len(Ns)):
        pyplot.plot(LDiff2, differenceData[i])
        pyplot.text(LDiff2[-1] + 0.25, differenceData[i][-1], "$N = " + str(Ns[i]) + "$", horizontalalignment='left', verticalalignment='center', fontsize=15)
    pyplot.show()

spectraFile = open("BigSpectraCollection.p", 'rb')
spectra = pickle.load(spectraFile)
spectraFile.close()

spectraFile2 = open("FractionalSprectra.p", 'rb')
spectraFrac = pickle.load(spectraFile2)
spectraFile2.close()
"""
N = 100
Ns = []
offDiags = []
keepGoing = True
while keepGoing:
    x, Zs = IQHEDiag.findEnergiesForRangeOfL(N, 4, 1, 0)
    Ns += [N]
    offDiags += [Zs[3][2][2]]
    N += 1
    if N == 350:
        keepGoing = False

(a, b), M = scipy.optimize.curve_fit(pow, Ns, offDiags)
print(a)
print(b)
pyplot.plot(Ns, offDiags)
pyplot.show()

#IQHEDiag.dumpRequest()

#spectrumCompareWithNoScatter(spectraFrac[(8, 6)], 8)
#scaleTest(spectra, f)
#peakAnalysis()
#differenceGraphPlotter()
#print(spectraFrac)
#for N in range(200, 301):


spectrum = FQHEDiag.findEnergiesForRangeOfL(5, 8, 3, 1, 0)
spectraFrac[(8, 6)] = spectrum
spectraFile2 = open("FractionalSprectra.p", 'wb')
pickle.dump(spectraFrac, spectraFile2)
spectraFile2.close()
"""


spectrum, Zs = IQHEDiag.findEnergiesForRangeOfL(50, 5, 1, 0)
spectrumCompareWithNoScatter(spectrum, 50)
