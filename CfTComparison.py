#CFT comparison
import math
import numpy
import numpy.linalg
import scipy
from usefulTools import generatePartitions
import matplotlib
import matplotlib.pyplot as pyplot
import IQHEDiag
import CFTMatrixElements

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
        return 1
    else:
        return doubleFactorial(n - 2)

def H22(state1, state2, maxOrder, N, m):
    answer = CFTMatrixElements.TOperatorMatrixElement((2,2), state1, state2)
    if state1 == state2:
        print(state1)
        print(answer)
        print("H22")
        print(" ")
    if maxOrder > 0:
        for i in range(1, maxOrder + 1):
            answer += 8*doubleFactorial(2*i + 1)*CFTMatrixElements.TOperatorMatrixElement((2,)*i, state1, state2)/doubleFactorial(2*i + 4)/(N*math.sqrt(m))**i
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
            answer += (alpha(i)*CFTMatrixElements.TOperatorMatrixElement((3,3) + (2,)*i, state1, state2) + beta(i)*CFTMatrixElements.TOperatorMatrixElement((2,)*(i+1), state1, state2))/(N*math.sqrt(m))**i
    return answer

def H(h_22, h_33, state1, state2, maxOrder, N, m):
    return h_22*H22(state1, state2, maxOrder, N, m) + h_33*H33(state1, state2, maxOrder, N, m)

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

def spectrumCompare(numericalSpectrum, maxOrder, N, m):
    LMax = numericalSpectrum[-1][0] + 1
    print(numericalSpectrum)
    E_0, h_22, h_33 = findParameters(numericalSpectrum)
    CFTSpectrum = calcSpetrum(LMax, E_0, h_22, h_33, maxOrder, N, m)
    L1 = [item[0] for item in numericalSpectrum]
    E1 = [item[1] for item in numericalSpectrum]
    L2 = [item[0] for item in CFTSpectrum]
    E2 = [item[1] for item in CFTSpectrum]
    pyplot.xlabel("Delta L")
    pyplot.ylabel("E/(e^2/epsilon0/magnetic length/(4*pi))")
    pyplot.plot(L1, E1, 'bo')
    pyplot.plot(L2, E2, 'rx')
    pyplot.show()

spectrumCompare(IQHEDiag.findEnergiesForRangeOfL(30, 6, 1, 0), 7, 30, 1)
