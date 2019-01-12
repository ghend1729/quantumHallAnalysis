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
import potentials
import haldanePotentials

mpmath.mp.dps = 25

matrixElementMemory = {}

useHaldane = True
V = potentials.exponentailRepulsion

def matrixElement(magneticLength, m1Prime, m2Prime, m1, m2):
    if (m1Prime, m2Prime, m1, m2) in matrixElementMemory:
        result = matrixElementMemory[(m1Prime, m2Prime, m1, m2)]
    elif (m1, m2, m1Prime, m2Prime) in matrixElementMemory:
        result = matrixElementMemory[(m1, m2, m1Prime, m2Prime)]
    else:
        if useHaldane:
            result = haldanePotentials.haldaneMatrixElement(m1Prime, m2Prime, m1, m2, magneticLength, V)
        else:
            result = matrixElementC(magneticLength, m1Prime, m2Prime, m1, m2)
        matrixElementMemory[(m1Prime, m2Prime, m1, m2)] = result
    return result

def NElectronMatrixElement(state1, state2, magneticLength):
    """
    Calculate the interaction matrix element between two slater determinent states.
    """
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
    return [state1Diff, state2Diff]
