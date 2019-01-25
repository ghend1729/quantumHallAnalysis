#CFT comparison
import math
import numpy
from numpy.linalg import inv
import numpy.linalg
import scipy
import scipy.integrate
from usefulTools import generatePartitions
import matplotlib
import matplotlib.pyplot as pyplot
import IQHEDiag
import CFTMatrixElements
import pickle

def alpha(n):
    if n == 0:
        return 1
    else:
        return (2*n + 3)*alpha(n-1)/(2*n)

def beta(n):
    if n == 0:
        return 0
    else:
        return (2*n + 1)*beta(n-1)/(2*n + 4) - 6*(2*n + 1)*alpha(n)/((n+1)*(n+2))

def doubleFactorial(n):
    if n == 1 or n == 2:
        return n
    else:
        return n*doubleFactorial(n - 2)

def H22(state1, state2, maxOrder, N, m):
    answer = CFTMatrixElements.TOperatorMatrixElement((2,2), state1, state2)
    if state1 == state2:
        print(state1)
        print(answer)
        print("H22")
        print(" ")
    if maxOrder > 0:
        for i in range(1, maxOrder + 1):
            answer += 8*doubleFactorial(2*i + 1)*CFTMatrixElements.TOperatorMatrixElement((2,)*(i+2), state1, state2)/doubleFactorial(2*i + 4)/((N*math.sqrt(m))**i)
    return answer

def H33(state1, state2, maxOrder, N, m):
    answer = CFTMatrixElements.TOperatorMatrixElement((3,3), state1, state2)
    if state1 == state2:
        print(state1)
        print(answer)
        print("H33")
        print(" ")
    if maxOrder > 0:
        for i in range(1, maxOrder + 1):
            answer += (alpha(i)*CFTMatrixElements.TOperatorMatrixElement((3,3) + (2,)*i, state1, state2) + beta(i)*CFTMatrixElements.TOperatorMatrixElement((2,)*(i+1), state1, state2))/((N*math.sqrt(m))**i)
    return answer

def H(h_22, h_33, state1, state2, maxOrder, N, m):
    return h_22*H22(state1, state2, maxOrder, N, m) + h_33*H33(state1, state2, maxOrder, N, m)

def zeroOrderEnergies(n, h_22, h_33, h_44=0, h_55=0):
    return 2*n*(n**2 - 1)*(h_33*(n**2 - 4) - h_22 - (n**2 - 4)*(n**2 - 9)*h_44 + (n**2 - 4)*(n**2 - 9)*(n**2 - 16)*h_55)

def stateEnergy(partition, h_22, h_33):
    return sum([zeroOrderEnergies(n, h_22, h_33) for n in partition])

def findParameters(spectrum):
    for i in spectrum:
        print(i)
    E_0 = spectrum[0][1]
    E_2 = spectrum[2][1]
    E_3 = spectrum[4][1]
    h_22 = (E_0 - E_2)/12
    h_33 = (E_3 - E_0 + 48*h_22)/240
    return E_0, h_22, h_33

def calcSpetrumZeroOrder(LMax, E_0, h_22, h_33):
    result = []
    for L in range(LMax):
        partitions = generatePartitions(L)
        result += [[L, E_0 + stateEnergy(p, h_22, h_33)] for p in partitions]
    return result

def calcSpetrum(LMax, E_0, h_22, h_33, maxOrder, N, m):
    result = []
    for L in range(LMax):
        partitions = generatePartitions(L)
        levelMatrix = numpy.array([[H(h_22, h_33, state1, state2, maxOrder, N, m) for state2 in partitions] for state1 in partitions])
        print(levelMatrix)
        energies = numpy.linalg.eigvals(levelMatrix)
        result += [[L, E_0 + E] for E in energies]
    return result

def coulombIntegrand(x, n):
    return math.cos(n*x)/math.sin(x/2)

def coulombIntegrandV2(x, n, N, m):
    return math.cos(n*x)/(2*math.sin(x/2) + 1/math.sqrt(2*N*m))

def coulombIntegrandV3(x, n, N, m):
    return math.cos(n*x)/(math.sqrt((2*math.sin(x/2))**2 + 1/2*N*m))

def xi(n, N, m):
    if n == 1:
        return 0
    else:
        y = scipy.integrate.quad(lambda x: coulombIntegrandV2(x, n, N, m), 0, math.pi)
        return n*y[0]

def stateEnergyCoulomb(p, N, m):
    return sum([xi(n, N, m) for n in p])

def spectrumCompare(numericalSpectrum, maxOrder, N, m):
    LMax = numericalSpectrum[-1][0] + 1
    print(numericalSpectrum)
    E_0, h_22, h_33 = findParameters(numericalSpectrum)
    print(h_22, h_33)
    CFTSpectrum = calcSpetrum(LMax, E_0, h_22, h_33, maxOrder, N, m)
    print(CFTSpectrum)
    L1 = [item[0] for item in numericalSpectrum]
    E1 = [item[1] for item in numericalSpectrum]
    U = 1.5*(max(E1) - min(E1))/(LMax - 1)
    print("U = " + str(U))
    E1 = [item[1] + U*item[0] for item in numericalSpectrum]
    L2 = [item[0] for item in CFTSpectrum]
    E2 = [item[1] + U*item[0] for item in CFTSpectrum]
    pyplot.xlabel("Angular momentum above ground state", fontsize = 24)
    pyplot.ylabel("Energy/w_0", fontsize = 24)
    pyplot.plot(L1, E1, 'bo')
    pyplot.hlines(E2, [i - 0.2 for i in L2], [i + 0.2 for i in L2])
    pyplot.show()

def calcSpetrumCoulomb(LMax, E_0, h_22, h_33, f, N, m):
    result = []
    for L in range(LMax):
        partitions = generatePartitions(L)
        result += [[L, E_0 + f*stateEnergyCoulomb(p, N, m) + stateEnergy(p, h_22, h_33)] for p in partitions]
    return result


def spectrumCompareCoulomb(numericalSpectrum, maxOrder, N, m):
    LMax = numericalSpectrum[-1][0] + 1
    print(numericalSpectrum)
    E_0 = numericalSpectrum[0][1]
    E_2 = numericalSpectrum[2][1]
    E_3 = numericalSpectrum[4][1]
    E_4 = numericalSpectrum[7][1]
    energyArray = [E_2 - E_0, E_3 - E_0, E_4 - E_0]
    r22 = [-2*n*(n**2 - 1) for n in range(2, 5)]
    r33 = [2*n*(n**2 - 4)*(n**2 - 1) for n in range(2, 5)]
    rXi = [xi(n, N, m) for n in range(2, 5)]
    fitterMatrix = numpy.transpose([r22, r33, rXi])
    coeficients = inv(fitterMatrix).dot(energyArray)
    h_22 = coeficients[0]
    h_33 = coeficients[1]
    f = coeficients[2]
    CFTSpectrum = calcSpetrumCoulomb(LMax, E_0, h_22, h_33, f, N, m)
    print(CFTSpectrum)
    L1 = [item[0] for item in numericalSpectrum]
    E1 = [item[1] for item in numericalSpectrum]
    U = 0*(max(E1) - min(E1))/(LMax - 1)
    print("U = " + str(U))
    E1 = [item[1] + U*item[0] for item in numericalSpectrum]
    L2 = [item[0] for item in CFTSpectrum]
    E2 = [item[1] + U*item[0] for item in CFTSpectrum]
    pyplot.xlabel("Angular momentum above ground state", fontsize = 24)
    pyplot.ylabel("Energy/(e^2/(4*pi*epsilon0*l_B)", fontsize = 24)
    pyplot.plot(L1, E1, 'bo')
    pyplot.hlines(E2, [i - 0.2 for i in L2], [i + 0.2 for i in L2])
    pyplot.show()

m = 1
N = 90
maxOrder = 1

spectraFile = open("Spectra.p", 'rb')
spectra = pickle.load(spectraFile)
spectraFile.close()

if N in spectra:
    spectrumCompareCoulomb(spectra[N], maxOrder, N, m)
else:
    spectrum = IQHEDiag.findEnergiesForRangeOfL(N, 8, 1, 0)
    spectra[N] = spectrum
    spectraFile2 = open("Spectra.p", 'wb')
    pickle.dump(spectra, spectraFile2)
    spectraFile2.close()
    spectrumCompare(spectra[N], maxOrder, N, m)
