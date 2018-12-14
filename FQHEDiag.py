#FQHE
import scipy
import numpy
import mpmath
import sympy
import itertools
import waveFunctionClasses
import DiagHamInterface
import matplotlib
import matplotlib.pyplot as pyplot
from usefulTools import generatePartitions as partitions
from usefulTools import signOfPermutation


def genStates(L, N, m, magneticLength):
    LaughlinState = [i*m for i in range(N)]
    if m % 2 == 1:
        fermion = True
    else:
        fermion = False
    levelPartitions = partitions(L)
    jackStates = []
    for p in levelPartitions:
        newJackState = LaughlinState[:]
        for i in range(len(p)):
            newJackState[N - i - 1] += p[len(p) - i - 1]
        jackStates.append(newJackState)
    decomposedStates = [DiagHamInterface.decomposeJackPolyState(p, fermion=fermion) for p in jackStates]
    jackBasis = [waveFunctionClasses.waveFunction(decompJack, magneticLength, fermion=fermion, convertToNormalisedBasis=True) for decompJack in decomposedStates]
    waveFunctions = waveFunctionClasses.gramSchmidt(jackBasis)
    return waveFunctions

def diagLevelL(L, N, m, magneticLength):
    states = genStates(L, N, m, magneticLength)
    numOfStates = len(states)

    halfMatrix = [[waveFunctionClasses.waveFuncMatrixElement(states[i], states[j]) for j in range(i+1)] for i in range(numOfStates)]
    #halfMatrix = [[longFormMatrixElement(magneticLength, states[i], states[j]) for j in range(i+1)] for i in range(numOfStates)]
    transposedHalfMatrix = [[halfMatrix[j][i] for j in range(i+1, numOfStates)] for i in range(numOfStates - 1)]
    print(transposedHalfMatrix)
    fullMatrix = [halfMatrix[i] + transposedHalfMatrix[i] for i in range(numOfStates - 1)]
    fullMatrix.append(halfMatrix[numOfStates-1])
    pertubationMatrix = mpmath.mp.matrix(fullMatrix)
    print(pertubationMatrix)
    energies = mpmath.mp.eigsy(pertubationMatrix, eigvals_only = True)
    return [float(mpmath.nstr(x)) for x in energies]

def findEnergiesForRangeOfL(LMax, N, m, magneticLength, alpha):
    finalList = []
    groundConfinementEnergy = alpha*N*m*(N-1)/2
    for L in range(LMax):
        finalList += [[L, E + groundConfinementEnergy + alpha*L] for E in diagLevelL(L, N, m, magneticLength)]
    return finalList

def plotEnergies(N, m, magneticLength, LMax, alpha):
    LEList = findEnergiesForRangeOfL(LMax, N, m, magneticLength, alpha)
    L = [item[0] for item in LEList]
    E = [item[1] + alpha*item[0] for item in LEList]
    pyplot.xlabel("Delta L")
    pyplot.ylabel("E/w_0")
    pyplot.plot(L, E, 'bo')
    pyplot.show()

plotEnergies(12, 3, 1, 5, 0.0596/12*3)
