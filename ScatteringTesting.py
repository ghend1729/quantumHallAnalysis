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

def f(n, g, h, a):
    return (h + g*n + a*(n**3))*(n-1)

def spectrumCompareWithNoScatter(numericalSpectrum, N):
    LMax = numericalSpectrum[-1][0] + 1
    print(numericalSpectrum)
    E_0 = numericalSpectrum[0][1]
    E = findDispersion(numericalSpectrum, LMax)
    CFTSpectrum = predictSpectrum(E, LMax, E_0)
    print(CFTSpectrum)
    L1 = [item[0] for item in numericalSpectrum if item[0] < 8]
    E1 = [item[1] for item in numericalSpectrum if item[0] < 8]
    U = 2.1*(max(E1) - min(E1))/(LMax - 1)
    print("U = " + str(U))
    E1 = [item[1] + U*item[0] for item in numericalSpectrum if item[0] < 8]
    L2 = [item[0] for item in CFTSpectrum if item[0] < 8]
    E2 = [item[1] + U*item[0] for item in CFTSpectrum if item[0] < 8]
    Ls = [L for L in E if L > 0]
    EAbs = [-E[L] for L in E if L > 0]
    LFit = [L for L in E if L > 0]
    EFit = [-E[L] for L in E if L > 0]
    g, h, a = scipy.optimize.curve_fit(f, LFit, EFit)[0]
    print(g, h, a)
    EPred = [f(L, g, h, a) for L in Ls]
    EDiff = [(EAbs[i] - EPred[i])/max(EAbs) for i in range(len(EAbs))]
    EDiff2 = [EAbs[i] - EAbs[i-1] for i in range(1, len(Ls))]
    LDiff2 = [Ls[i] for i in range(1, len(Ls))]
    findPeak(LDiff2, EDiff2)

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

    ax = pyplot.subplot(121)
    ax.tick_params(labelsize = 15)
    pyplot.title("$N = " + str(N) + "$ $m=3$ Free Boson Test", fontsize = 22)
    pyplot.xlabel("$\Delta L$", fontsize = 20)
    pyplot.ylabel("$\Delta E/(q^2/4\pi\epsilon_0l_B)$", fontsize = 20)
    pyplot.plot(L1, E1, 'ko', label = "$\delta\hat{H}$")
    pyplot.hlines(E2, [i - 0.3 for i in L2], [i + 0.3 for i in L2], label = "$\hat{H}$")
    pyplot.text(0.98, 0.02, "$U_0 = " + str(round(U/2, 4)) + "q^2/4\pi\epsilon_0l_B^3$", horizontalalignment='right', fontsize = 20, transform=ax.transAxes)
    pyplot.legend(fontsize = 22)

    ax2 = pyplot.subplot(122)
    ax2.tick_params(labelsize=15)
    pyplot.title("$m = 1$ Error vs N", fontsize=25)
    pyplot.xlabel("N", fontsize=22)
    pyplot.ylabel("ERROR(%)", fontsize = 22)
    pyplot.plot(Ns, errorsList, 'ko', label = "ACTUAL ERROR")
    pyplot.plot(NsPred, errorListPred, label = "ERROR = $bN^a$ FIT")
    pyplot.text(0.98, 0.75, "a = " + str(round(a, 3)) + " ± {0:.3f}".format(aError) + "\n" + "b = " + str(round(b, 2)) + " ± {0:.1f}".format(bError), horizontalalignment='right', verticalalignment='center', fontsize=22, transform=ax2.transAxes)
    pyplot.legend(fontsize = 22)
    pyplot.show()

def pow(x, a, b):
    return b*x**a

def scaleTest(spectra, f):
    Ns = []
    hs = []
    gs = []
    As = []
    asScaled = []
    for N in range(150, 340):
        spectrum = spectra[(N, 11)]
        E = findDispersion(spectrum, 11)
        Ls = [L for L in E if L > 0]
        EAbs = [-E[L] for L in E if L > 0]
        g, h, a = scipy.optimize.curve_fit(f, Ls, EAbs)[0]
        Ns.append(N)
        hs.append(h)
        gs.append(g)
        As.append(a)
        if N == 250:
            ax = pyplot.subplot(1,3,1)
            pyplot.plot(Ls, EAbs, 'ko')
            pyplot.plot(Ls, [f(L, g, h, a) for L in Ls])

    a, b = scipy.optimize.curve_fit(pow, Ns, As)[0]
    c, d = scipy.optimize.curve_fit(pow, Ns, gs)[0]
    e, f = scipy.optimize.curve_fit(pow, Ns, hs)[0]
    print(a, b)
    print(c, d)
    print(e, f)

    pyplot.subplot(1,3,2)
    pyplot.xlabel("N", fontsize = 20)
    pyplot.ylabel("g")
    pyplot.plot(Ns, gs)

    pyplot.subplot(1,3,3)
    pyplot.xlabel("N", fontsize = 20)
    pyplot.ylabel("h")
    pyplot.plot(Ns, hs)

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

spectraFile2 = open("FractionalSprectra.p", 'rb')
spectraFrac = pickle.load(spectraFile2)
spectraFile2.close()
"""
N = 200
keepGoing = True
while keepGoing:
    x = IQHEDiag.findEnergiesForRangeOfL(N, 11, 1, 0)
    y = spectra[(N,11)]
    for i in range(len(x)):
        if not (x[i][1] - y[i][1] == 0.0):
            keepGoing = False
            print("Error!!!")
    N += 1
    if N == 340:
        keepGoing = False

IQHEDiag.dumpRequest()
"""
spectrumCompareWithNoScatter(spectraFrac[(8, 6)], 8)
#scaleTest(spectra, f)
#peakAnalysis()
#print(spectraFrac)
#for N in range(200, 301):

"""
spectrum = FQHEDiag.findEnergiesForRangeOfL(6, 8, 3, 1, 0)
spectraFrac[(8, 6)] = spectrum
spectraFile2 = open("FractionalSprectra.p", 'wb')
pickle.dump(spectraFrac, spectraFile2)
spectraFile2.close()
"""
