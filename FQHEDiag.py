#FQHE
import scipy
import numpy
import mpmath
import sympy
import itertools
import waveFunctionClasses
import symPoly
import matplotlib
import matplotlib.pyplot as pyplot
from usefulTools import generatePartitions as partitions
from usefulTools import signOfPermutation

mpmath.mp.dps = 20

#number of particles
N = 5

#1/v
m = 3

#max change in angular momentum
LMax = 4

alpha = 0

magneticLength = 1

isFermion = True
#first createB lauglin functioin squared in the symetric power sum basis
baseWaveFunctionPoly = symPoly.genralSymPoly([])
newPartition = []
for p in itertools.permutations(range(N)):
    newPartition = [p[i] + i for i in range(N)]
    if p[0] == 0:
        baseWaveFunctionPoly = baseWaveFunctionPoly + symPoly.symetricPowerSumPoly(newPartition, N*signOfPermutation(p))
    else:
        baseWaveFunctionPoly = baseWaveFunctionPoly + symPoly.symetricPowerSumPoly(newPartition, signOfPermutation(p))

if False:
    #boson state1
    isFermion = False
    baseWaveFunctionPoly = baseWaveFunctionPoly**(m//2)
else:
    baseWaveFunctionPoly = baseWaveFunctionPoly**((m-1)//2)


def genStates(L):
    statePolys = [baseWaveFunctionPoly*symPoly.symetricPowerSumPoly(p, 1) for p in partitions(L)]
    if isFermion:
        waveFunctions = [waveFunctionClasses.waveFunctionFermion(p, N, magneticLength) for p in statePolys]
    return waveFunctions

def diagLevelL(L):
    states = genStates(L)
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

def findEnergiesForRangeOfL():
    finalList = []
    groundConfinementEnergy = alpha*N*m*(N-1)/2
    for L in range(LMax):
        finalList += [[L, E + groundConfinementEnergy + alpha*L] for E in diagLevelL(L)]
    return finalList

def plotEnergies():
    LEList = findEnergiesForRangeOfL()
    L = [item[0] for item in LEList]
    E = [item[1] for item in LEList]
    pyplot.xlabel("Delta L")
    pyplot.ylabel("E/(e^2/epsilon0/magnetic length/(4*pi))")
    pyplot.plot(L, E, 'bo')
    pyplot.show()

plotEnergies()
