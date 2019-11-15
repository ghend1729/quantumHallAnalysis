#Diagonalise the pertubation
"""
This entire module is for calculating matrix elements <A|V|B> when |A> and |B> are
slater determinents of lowest Landau level single particle state.

This code could be cleaned up perhaps but it works. As far as I can tell.
"""
import math
import numpy
import numpy.linalg
import scipy
import scipy.linalg
import mpmath
import pickle
import itertools
from sympy.combinatorics.permutations import Permutation
import matplotlib
import matplotlib.pyplot as pyplot
from CoulombMatrixFunctions import *
import potentials
import haldanePotentials

mpmath.mp.dps = 20
#Set mpmath percision
useHaldane = False
#indicate to wether to use Haldane pseudo-potentials. If this is false the coulomb matrix elements
#matrixElementC are used.
V = potentials.exponentailRepulsion

#All previously used matrix elements of the coulomb interaction have been stored to speed things up.
#The next two functions loads and updates this file. dumpMatrixElements is triggered by a function IQHEDiag.
if not useHaldane:
    matrixElementFile = open("coulombMatrixElements.p", 'rb')
    matrixElementMemory = pickle.load(matrixElementFile)
    matrixElementFile.close()
else:
    matrixElementMemory = {}

def dumpMatrixElements():
    """
    Use this to save current matrix elements. This can also be triggered from
    IQHEDiag module by using dumpRequest().
    """
    matrixElementsFile = open("coulombMatrixElements.p", 'wb')
    pickle.dump(matrixElementMemory, matrixElementsFile)
    matrixElementsFile.close()

def matrixElement(magneticLength, m1Prime, m2Prime, m1, m2):
    """
    Calculates <m1Prime|<m2Prime|V|m1>|m2>. Depending on wether V is the coulomb interaction
    or any other rotationally symmetric potential, it will use the haldane pseudo-potentails or use
    analytic coulomb matrix elements.

    It also checks if these have been calculated before and collects and strores
    past calculations.
    """
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
    This uses standard formulas that depends on wether state1 or state2 differ by 0, 1, 2 single particle states.
    If state1 and state2 differ by more than one single particle state then the matrix element is 0.
    The sign factors seen below are due to reordering of slater determinent states. We essentailly need
    the single particle states that state1 and state2 differ by are in the same 'slots' of the slater determinent.
    """
    exchangeFactor = -1
    x = 0
    #x will store the final answer.
    diff = findDifferentElements(state1, state2)
    if len(diff[0]) == 0:
        for i in range(len(state1)):
            for j in range(i+1, len(state1)):
                x += twoElectronMatrixElement(state1[i], state1[j], state2[i], state2[j], magneticLength, exchangeFactor)
    elif len(diff[0]) == 1:
        state1Diff = diff[0][0]
        state2Diff = diff[1][0]
        N = len(state1)
        newState1 = state1
        newState2 = state2
        del newState1[state1Diff]
        del newState2[state2Diff]

        sighn1 = (exchangeFactor)**(N-1 - state1Diff)
        sighn2 = (exchangeFactor)**(N-1 - state2Diff)

        x = sighn1*sighn2*sum([matrixElement(magneticLength, state2[i], state1[state1Diff], state2[i], state2[state2Diff]) + exchangeFactor*matrixElement(magneticLength, state1[state1Diff], state2[i], state2[i], state2[state2Diff]) for i in range(N)])
    elif len(diff[0]) == 2:
        N = len(state1)
        i, j = diff[0][0], diff[0][1]
        k, l = diff[1][0], diff[1][1]
        sighn1 = (exchangeFactor)**(2*N - 3 - i - j)
        sighn2 = (exchangeFactor)**(2*N - 3 - k - l)
        x = sighn1*sighn2*(matrixElement(magneticLength,state1[i],state1[j],state2[k],state2[l]) + exchangeFactor*matrixElement(magneticLength,state1[j],state1[i],state2[k],state2[l]))
    else:
        x = 0
    return x

def twoElectronMatrixElement(i,j,k,l, magneticLength, exchangeFactor):
    """
    Calculates matrix element between two correctly anti-symmerised two electrons states.
    """
    y = matrixElement(magneticLength, i, j, k, l) + exchangeFactor*matrixElement(magneticLength, i, j, l, k)
    return y

def findDifferentElements(state1, state2):
    """
    Finds the single particle states taht state1 and state2 differ by.
    """
    state1Diff = [i for i in range(len(state1)) if not state1[i] in state2]
    state2Diff = [i for i in range(len(state1)) if not state2[i] in state1]
    return [state1Diff, state2Diff]
