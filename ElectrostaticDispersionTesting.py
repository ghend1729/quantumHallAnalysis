#Electrostatic Dispersion Testing
import numpy
import scipy
import scipy.special as special
import scipy.integrate as integrate
import math
import mpmath
import IQHEDiag
import FQHEDiag
import matplotlib
import matplotlib.pyplot as pyplot
import usefulTools
import scipy.optimize
import pickle
import ScatteringTesting

mpmath.mp.dps = 10

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

def besselIntegrand(r, r0, k):
    return -k*special.j0(r*k)*special.j1(r0*k)

def integrateTest(theLimit, r, r0):
    return integrate.quad(lambda x: besselIntegrand(r, r0, x), 0, theLimit, limit = 1000)

def besselIntegrandReversed(r, k, r0, N):
    return r*special.j0(k*r)*chargeDistro(r*r0, N)

def kIntegrand(r0, k, N):
    return -2*math.pi*k*special.j1(k)*integrate.quad(lambda r: besselIntegrandReversed(r, k, r0, N), 0, numpy.inf, limit = 1000)[0]

def integrateBessel(r, r0):
    """
    compute the bessel integrals in the case special averaging limits.
    """
    delta = abs(r - r0)
    if delta == 0:
        delta = 2*math.pi/(r + r0)
    else:
        delta = 2*math.pi/delta
    chunks = 7
    return sum([integrateTest(2000 + i*delta/chunks, r, r0)[0] for i in range(chunks)])/chunks

def chargeDistro(r, N):
    R = math.sqrt(2*N)
    if r < R:
        return 1/(2*math.pi)
    else:
        return numpy.exp(-(r-R)**2)/(2*math.pi)

def chargeDistroExact(r, N):
    answer = sum([mpmath.power((r**2)/2, m)/mpmath.factorial(m) for m in range(N-1)])*mpmath.exp(-(r**2)/2)/(2*math.pi)
    return float(mpmath.nstr(answer, n=10))


def plotPotential():
    x = [1/i for i in range(10000, 10040)]
    y = [differentiator(unnormalisedPotential, 1, r) for r in x]
    pyplot.xlabel("r", fontsize = 18)
    pyplot.ylabel("V(r)", fontsize = 18)
    pyplot.plot(x, y)
    pyplot.show()

def differentiator2(r0, R):
    dx = 0.01
    x = [i*dx for i in range(1000)]
    y = [2*math.pi*integrateBessel(r, 1)*r*chargeDistro(r*r0, R) for r in x]
    return dx*sum(y)

def differentiator3(r0, N):
    return integrate.quad(lambda k: kIntegrand(r0, k, N), 0, numpy.inf, limit = 1000)[0]

def spectrumFitter(spectrum, N):
    LMax = spectrum[-1][0] + 1
    E = ScatteringTesting.findDispersion(spectrum, LMax)
    EAbs = [-E[L] for L in E if L > 0]
    Ls = [L for L in E]
    R = math.sqrt(2*(N))
    h = -differentiator3(R, N)/R
    p = (EAbs[13] - EAbs[0])/13
    print(N, p/h)
    """
    predE = [L*h for L in Ls]
    predE2 = [p[0]*L for L in Ls]
    pyplot.xlabel("n", fontsize = 18)
    pyplot.ylabel("|E(n)|", fontsize = 18)
    pyplot.plot(Ls, EAbs, 'ko')
    pyplot.plot(Ls, predE)
    pyplot.plot(Ls, predE2)
    pyplot.show()
    """
for N in [30, 60, 90]:
    spectrumFitter(ScatteringTesting.spectra[(N, 15)], N)


#plotPotential()
