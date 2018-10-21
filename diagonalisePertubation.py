#Diagonalise the pertubation
import math
import numpy
import numpy.linalg
import scipy
import matplotlib
import matplotlib.pyplot as pyplot
from CoulombMatrixFunctions import *

matrixElementMemory = {}

def NElectronMatrixElement(state1, state2, magneticLength):
    x = 0
    for i in range(len(state1)):
        for j in range(i+1, len(state1)):
            if includeInPotentialEnergySum(state1, state2, i, j, len(state1)):
                if (state1[i], state1[j], state2[i], state2[j]) in matrixElementMemory:
                    y = matrixElementMemory[(state1[i], state1[j], state2[i], state2[j])]
                else:
                    y = matrixElement(magneticLength, state1[i], state1[j], state2[i], state2[j]) - matrixElement(magneticLength, state1[i], state1[j], state2[j], state2[i])
                    matrixElementMemory[(state1[i], state1[j], state2[i], state2[j])] = y
                x += y
    return x

def includeInPotentialEnergySum(state1, state2, i, j, N):
    return all([state1[k] == state2[k] for k in range(N) if not (k==i or k==j)])

def generatePartitions(L):
    #source: https://stackoverflow.com/questions/10035752/elegant-python-code-for-integer-partitioning
    answer = set()
    answer.add((L, ))
    for x in range(1, L):
        for y in generatePartitions(L - x):
            answer.add(tuple(sorted((x, ) + y)))
    return answer

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
    print(halfMatrix)
    transposedHalfMatrix = [[halfMatrix[j][i] for j in range(i+1, numOfStates)] for i in range(numOfStates - 1)]
    print(transposedHalfMatrix)
    fullMatrix = [halfMatrix[i] + transposedHalfMatrix[i] for i in range(numOfStates - 1)]
    fullMatrix.append(halfMatrix[numOfStates-1])
    print(fullMatrix)
    pertubationMatrix = numpy.array(fullMatrix)
    print(pertubationMatrix)
    return numpy.linalg.eigvals(pertubationMatrix)/magneticLength

def findEnergiesForRangeOfL(N, LMax, magneticLength, alpha):
    finalList = []
    for L in range(LMax):
        finalList += [[L, E + alpha*L] for E in diagonaliseLLevel(L, N, magneticLength)]
    L = [item[0] for item in finalList]
    E = [item[1] for item in finalList]

    return L, E

def plotEnergies(N, LMax, magneticLength, alpha):
    L, E = findEnergiesForRangeOfL(N, LMax, magneticLength, alpha)
    pyplot.plot(L, E, 'bo')
    pyplot.show()

plotEnergies(80, 12, 1, 1/18)
