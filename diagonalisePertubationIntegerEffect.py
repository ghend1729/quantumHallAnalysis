#Diagonalise the pertubation
import math
import numpy
import numpy.linalg
import scipy
import scipy.linalg
import mpmath
import itertools
from sympy.combinatorics.permutations import Permutation
import matplotlib
import matplotlib.pyplot as pyplot
from CoulombMatrixFunctions import *
from usefulTools import generatePartitions

mpmath.mp.dps = 50

NElectronMatrixElementMemory = {}

matrixElementMemory = {}

def NElectronMatrixElement(state1, state2, magneticLength):
    if (state1, state2, magneticLength) in NElectronMatrixElementMemory:
        x = NElectronMatrixElementMemory[(state1, state2, magneticLength)]
    elif (state2, state1, magneticLength):
        x = NElectronMatrixElementMemory[(state2, state1, magneticLength)]
    else:
        x = 0
        diff = findDifferentElements(state1, state2)
        if len(diff[0]) == 0:
            for i in range(len(state1)):
                for j in range(i+1, len(state1)):
                    x += twoElectronMatrixElement(state1[i], state1[j], state2[i], state2[j], magneticLength)
        elif len(diff[0]) == 1:
            state1Diff = diff[0][0]
            state2Diff = diff[1][0]
            N = len(state1)
            newState1 = state1
            newState2 = state2
            del newState1[state1Diff]
            del newState2[state2Diff]

            sighn1 = (-1)**(N-1 - state1Diff)
            sighn2 = (-1)**(N-1 - state2Diff)

            x = sighn1*sighn2*sum([matrixElement(magneticLength, state2[i], state1[state1Diff], state2[i], state2[state2Diff]) - matrixElement(magneticLength, state1[state1Diff], state2[i], state2[i], state2[state2Diff]) for i in range(N)])
        elif len(diff[0]) == 2:
            N = len(state1)
            i, j = diff[0][0], diff[0][1]
            k, l = diff[1][0], diff[1][1]
            sighn1 = (-1)**(2*N - 3 - i - j)
            sighn2 = (-1)**(2*N - 3 - k - l)
            x = sighn1*sighn2*(matrixElement(magneticLength,state1[i],state1[j],state2[k],state2[l]) - matrixElement(magneticLength,state1[j],state1[i],state2[k],state2[l]))
        else:
            x = 0
        NElectronMatrixElementMemory[(state1, state2, magneticLength)] = x

    return x

def perm_parity(lst):
    return Permutation(lst).signature()

def permuteList(perm, listInput):
    return [listInput[i] for i in perm]

def potentialij(i, j, state1, state2, magneticLength):
    if all([state1[k] == state2[k] for k in range(len(state1)) if not (k==i or k==j)]):
        return matrixElement(magneticLength, state1[i], state1[j], state2[i], state2[j])
    else:
        return 0

def longFormMatrixElement(magneticLength, state1, state2):
    perms = itertools.permutations(range(len(state1)))
    x = 0
    N = len(state1)
    for p1 in perms:
        for p2 in perms:
            sighn = perm_parity(p1)*perm_parity(p2)
            for i in range(N):
                for j in range(i+1, N):
                    x += sighn*potentialij(i, j, permuteList(p1, state1), permuteList(p2, state2), magneticLength)
    return x

def twoElectronMatrixElement(i,j,k,l, magneticLength):
    y = 0
    if (i, j, k, l) in matrixElementMemory:
        y = matrixElementMemory[(i, j, k, l)]
    else:
        y = matrixElement(magneticLength, i, j, k, l) - matrixElement(magneticLength, i, j, l, k)
        matrixElementMemory[(i, j, k, l)] = y
    return y


def findDifferentElements(state1, state2):
    state1Diff = [i for i in range(len(state1)) if not state1[i] in state2]
    state2Diff = [i for i in range(len(state1)) if not state2[i] in state1]
    print(state1[70:])
    print(state2[70:])
    print(state1Diff)
    print(state2Diff)
    return [state1Diff, state2Diff]

def generateStates(L, N):
    partitions = [item for item in generatePartitions(L) if len(item) <= N]
    states = []
    for x in partitions:
        tempState = [i for i in range(N)]
        y = len(x)
        for i in range(y):
            tempState[N-1-i] = tempState[N-1-i] + x[y-1-i]
        states.append(tempState)
    return states

def diagonaliseLLevel(L,N, magneticLength):
    states = generateStates(L,N)
    numOfStates = len(states)

    halfMatrix = [[NElectronMatrixElement(states[i], states[j], magneticLength) for j in range(i+1)] for i in range(numOfStates)]
    #halfMatrix = [[longFormMatrixElement(magneticLength, states[i], states[j]) for j in range(i+1)] for i in range(numOfStates)]
    transposedHalfMatrix = [[halfMatrix[j][i] for j in range(i+1, numOfStates)] for i in range(numOfStates - 1)]
    print(transposedHalfMatrix)
    fullMatrix = [halfMatrix[i] + transposedHalfMatrix[i] for i in range(numOfStates - 1)]
    fullMatrix.append(halfMatrix[numOfStates-1])
    pertubationMatrix = mpmath.mp.matrix(fullMatrix)
    print(pertubationMatrix)
    energies = mpmath.mp.eigsy(pertubationMatrix, eigvals_only = True)
    return [float(mpmath.nstr(x)) for x in energies]


def findEnergiesForRangeOfL(N, LMax, magneticLength, alpha):
    finalList = []
    groundConfinementEnergy = alpha*N*(N-1)/2
    for L in range(LMax):
        finalList += [[L, E + groundConfinementEnergy + alpha*L] for E in diagonaliseLLevel(L, N, magneticLength)]
    return finalList

def plotEnergies(N, LMax, magneticLength, alpha):
    LEList = findEnergiesForRangeOfL(N, LMax, magneticLength, alpha)
    L = [item[0] for item in LEList]
    E = [item[1] for item in LEList]
    pyplot.xlabel("Delta L")
    pyplot.ylabel("E/(e^2/epsilon0/magnetic length/(4*pi))")
    pyplot.plot(L, E, 'bo')
    pyplot.show()

#plotEnergies(90, 8, 1, 1/18)
