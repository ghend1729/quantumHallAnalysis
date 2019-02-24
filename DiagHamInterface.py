#DiagHam Python Interface
"""
This is a very basic interface for the program to call DiagHam in order to obtain slater
or monomial decomposition of edge states.

This essentailly asks for slater decomp of jack polynomials S_lambda where lambda is a
partition that labels the jack polynomial.

When inputting partitions the number of elements should match the number of particles
so if the number of non-zero entries does not match zeros should be added until it does.
"""

import numpy
import os
import ast

def convertPartionForDiagHam(partition):
    """
    Converts the format of a partition to the one required for DiagHam to
    understand.
    """
    LzMax = max(partition)
    OccupationNumList = [0 for i in range(LzMax + 1)]
    for i in range(LzMax + 1):
        for x in partition:
            if i == x:
                OccupationNumList[i] += 1

    resultString = ""
    for i in range(LzMax):
        resultString += str(OccupationNumList[i]) + " "
    resultString += str(OccupationNumList[-1])
    return resultString

def makeStateReferenceFile(partition):
    """
    Makes input file for DiagHam when we wish a slater decomp of
    a jack polynomial given by lambda = partition.
    """
    NbrParticles = len(partition)
    LzMax = max(partition)
    ReferenceState = convertPartionForDiagHam(partition)
    f = open("edgeState.dat", "w")
    f.write("NbrParticles=" + str(NbrParticles) + "\n")
    f.write("LzMax=" + str(LzMax) + "\n")
    f.write("ReferenceState=" + ReferenceState)
    f.close()

def getNBodyBasisDecomposition(partition, fermion = True):
    """
    Promps DiagHam to slater decomp jack polynomial with lambda = partition
    and output to a test file.
    """
    makeStateReferenceFile(partition)
    if fermion:
        os.system("../DiagHam/build//FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -a -2 -t edgeStateDecomp.txt --fermion --reference-file edgeState.dat")
    else:
        os.system("../DiagHam/build//FQHE/src/Programs/FQHEOnSphere/FQHESphereJackGenerator -a -2 -t edgeStateDecomp.txt --reference-file edgeState.dat")

def readInState():
    """
    Takes the output file from DiagHam and converts the decomp into a list containing each
    identifier and coeficient for each slater state.
    """
    f = open("edgeStateDecomp.txt", "r")
    fileLines = f.readlines()
    finalState = []
    for l in fileLines:
        splitLine = l.split(" ")
        coeficient = float(splitLine[0])
        nBodyState = sorted(ast.literal_eval(splitLine[1]))
        finalState.append([coeficient, nBodyState])
    return finalState

def decomposeJackPolyState(partition, fermion=True):
    """
    Uses diagham to slater decompose a jackpolynomial and returns this as a list containing each
    identifier and coeficient for each slater state.
    """
    getNBodyBasisDecomposition(partition, fermion)
    state = readInState()
    os.system("rm edgeState.dat")
    os.system("rm edgeStateDecomp.txt")
    particleNum = len(partition)
    for i in range(len(state)):
        if not len(state[i][1]) == particleNum:
            for j in range(state[i][1] - particleNum):
                state[i][1] = [0] + state[i][1]
        state[i][1] = tuple(state[i][1])
    return state
