#Wave Function Classes
import math
import scipy
import numpy
import diagonalisePertubationIntegerEffect

def singleParticleNorm(n, magneticLength):
    return 1/(math.sqrt(math.pi*math.factorial(n))*magneticLength**(n+1))

def slaterNorm(state, magneticLength):
    n = len(state)
    answer = 1/math.sqrt(math.factorial(n))
    for m in state:
        answer = answer*singleParticleNorm(m, magneticLength)
    return answer

def convertPartitionToState(partition, n):
    baseState = [i for i in range(n)]
    for i in range(len(partition)):
        baseState[n - len(partition) + i] += partition[i]
    return baseState

class waveFunctionFermion(poly, n, magneticLength):
    def __init__(self, poly, n, magneticLength):
        schurConverted = poly.convertToShurBasis()
        self.states = [[s.coeficient, convertPartitionToState(s.partition, n)] for s in schurConverted.schurs]
        self.n = n
        self.magneticLength = magneticLength
        self.normalise()

    def normalise(self):
        for i in range(len(self.states)):
            self.states[i][0] = self.states[i][0]/slaterNorm(self.states[i][1], magneticLength)
        normConst = math.sqrt(sum([(x[0])**2 for x in self.states]))
        for i in range(len(self.states)):
            self.states[i][0] = self.states[i][0]/normConst

def waveFuncMatrixElement(state1, state2):
    answer = 0
    for s1 in state1.states:
        for s2 in state2.states:
            answer += diagonalisePertubationIntegerEffect.NElectronMatrixElement(s1[1], s2[1], state1.magneticLength)*s1[0]*s2[0]
    return answer
