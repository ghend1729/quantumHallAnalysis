#Integer Quantum Hall Effect
"""
This module is for calculating the edge spectrum of the integer quantum hall effect.
"""
import math
import numpy
import numpy.linalg
import scipy
import scipy.linalg
import mpmath
import matplotlib
import matplotlib.pyplot as pyplot
import NBodyBasisMatrixElementCalc
from usefulTools import generatePartitions

def dumpRequest():
    """
    Save the matrix element memeroy.
    """
    NBodyBasisMatrixElementCalc.dumpMatrixElements()

def generateStates(L, N):
    """
    Generate a slater basis of states for angular momentum level L above ground.
    """
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
    """
    Calculates the pertubation matrix for angular momentum L above
    the Laughlin state and then diagonalises it to find the first order perturbation
    to the energy levels (when we have N particles).
    Returns these energies as a list.
    """
    states = generateStates(L,N)
    numOfStates = len(states)

    halfMatrix = [[NBodyBasisMatrixElementCalc.NElectronMatrixElement(states[i], states[j], magneticLength) for j in range(i+1)] for i in range(numOfStates)]
    #halfMatrix = [[longFormMatrixElement(magneticLength, states[i], states[j]) for j in range(i+1)] for i in range(numOfStates)]
    transposedHalfMatrix = [[halfMatrix[j][i] for j in range(i+1, numOfStates)] for i in range(numOfStates - 1)]
    print(transposedHalfMatrix)
    fullMatrix = [halfMatrix[i] + transposedHalfMatrix[i] for i in range(numOfStates - 1)]
    fullMatrix.append(halfMatrix[numOfStates-1])
    pertubationMatrix = mpmath.mp.matrix(fullMatrix)
    print(pertubationMatrix)
    print("Diagonalising L = " + str(L) + " level with N = " + str(N))
    energies = mpmath.mp.eigsy(pertubationMatrix, eigvals_only = True, overwrite_a = True)
    return [float(mpmath.nstr(x, n=20)) for x in energies]


def findEnergiesForRangeOfL(N, LMax, magneticLength, alpha):
    """
    Finds the first order pertubation energies for angular momentum up to LMax - 1 above
    Laughlin state. This will return a list of two element lists.
    The two element lists are in the form:
    [angular momentum above Laughlin state, energy above Laughlin state]
    """
    finalList = []
    groundConfinementEnergy = alpha*N*(N-1)/2
    for L in range(LMax):
        finalList += [[L, E + groundConfinementEnergy + alpha*L] for E in diagonaliseLLevel(L, N, magneticLength)]
    return finalList

def plotEnergies(N, LMax, magneticLength, U0):
    alpha = U0/N
    LEList = findEnergiesForRangeOfL(N, LMax, magneticLength, alpha)
    L = [item[0] for item in LEList]
    E = [item[1] for item in LEList]
    pyplot.xlabel("Delta L")
    pyplot.ylabel("E/(e^2/epsilon0/magnetic length/(4*pi))")
    pyplot.plot(L, E, 'bo')
    pyplot.show()
