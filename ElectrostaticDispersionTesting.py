#Electrostatic Dispersion Testing
import numpy
import scipy
import scipy.special as special
import scipy.integrate as integrate
import math
import IQHEDiag
import FQHEDiag
import matplotlib
import matplotlib.pyplot as pyplot
import usefulTools
import scipy.optimize
import pickle
import ScatteringTesting

def unnormalisedPotential(x):
    return integrate.quad(lambda q: (special.j0(q*x)/q)*special.j1(q), 0, numpy.inf, limit = 1000)[0]

def flippedPotential(x):
    return unnormalisedPotential(1) - unnormalisedPotential(x)

def xConverter(n, N):
    return math.sqrt((N - 1 + n)/(N - 1))

def derivativeOfPotential(x):
    return integrate.quad(lambda q: special.jvp(0,q*x)*special.j1(q), 0, numpy.inf, limit = 1000)[0]

def differentiator(f, x, dx):
    return (f(x+dx) - f(x))/dx

def plotPotential():
    x = [1/i for i in range(10000, 10040)]
    y = [differentiator(unnormalisedPotential, 1, r) for r in x]
    pyplot.xlabel("r", fontsize = 18)
    pyplot.ylabel("V(r)", fontsize = 18)
    pyplot.plot(x, y)
    pyplot.show()

def spectrumFitter(spectrum, N):
    LMax = spectrum[-1][0] + 1
    E = ScatteringTesting.findDispersion(spectrum, LMax)
    EAbs = [-E[L] for L in E]
    Ls = [L for L in E]
    h = -differentiator(unnormalisedPotential, 1, 1/1000000000)/math.sqrt(2*(N))
    p = numpy.polyfit(Ls[5:], EAbs[5:], 1)
    print(h/p[0])
    predE = [L*h for L in Ls]
    predE2 = [p[0]*L for L in Ls]
    pyplot.xlabel("n", fontsize = 18)
    pyplot.ylabel("|E(n)|", fontsize = 18)
    pyplot.plot(Ls, EAbs, 'ko')
    pyplot.plot(Ls, predE)
    pyplot.plot(Ls, predE2)
    pyplot.show()

spectrumFitter(ScatteringTesting.spectra[(339, 11)], 339)
#plotPotential()
