#Wave Function Classes
import math
import scipy
import numpy
import diagonalisePertubationIntegerEffect
import copy

def singleParticleNorm(n, magneticLength):
    return 1/(math.sqrt(math.pi*math.factorial(n)*2**(n+1))*magneticLength**(n+1))

def NBodyNorm(state, magneticLength):
    n = len(state)
    occupiedStates = set(state)
    answer = 1/math.sqrt(math.factorial(n))
    for m in occupiedStates:
        occupancyNum = state.count(n)
        answer = answer*(1/math.sqrt(math.factorial(occupancyNum)))*((singleParticleNorm(m, magneticLength))**occupancyNum)
    return answer

def convertPartitionToState(partition, n):
    baseState = [i for i in range(n)]
    for i in range(len(partition)):
        baseState[n - len(partition) + i] += partition[i]
    return baseState

class waveFunction:
    def __init__(self, statesDecomp, magneticLength, fermion=True, convertToNormalisedBasis=False, normalise=False):
        self.states = copy.deepcopy(statesDecomp)
        self.magneticLength = magneticLength
        self.fermion = fermion
        if convertToNormalisedBasis:
            for i in range(len(self.states)):
                self.states[i][0] = self.states[i][0]/NBodyNorm(self.states[i][1], self.magneticLength)

    def __add__(self, otherWaveFunction):
        answerList = []
        for state in otherWaveFunction.states:
            correspondingState = next((s for s in self.states if s[1] == state[1]), [0, state[1]])
            answerList.append([correspondingState[0] + state[0], state[1]])
        NBodyKetsInOtherWave = [s[1] for s in otherWaveFunction.states]
        leftOverList = [s for s in self.states if not (s[1] in NBodyKetsInOtherWave)]
        answerList = answerList + leftOverList
        return waveFunction(answerList, self.magneticLength, fermion=self.fermion)

    def __mul__(self, otherNum):
        answerList = [[s[0]*otherNum, s[1]] for s in self.states]
        return waveFunction(answerList, self.magneticLength, fermion=self.fermion)

    def __sub__(self, otherWaveFunction):
        return self + otherWaveFunction*(-1)

    def __or__(self, otherWaveFunction):
        answer = 0
        for state in otherWaveFunction.states:
            correspondingComponent = next((s[0] for s in self.states if s[1] == state[1]), 0)
            answer += state[0]*correspondingComponent
        return answer

    def normalise(self):
        sizeOfState = (self | self)
        self = self*(1/math.sqrt(sizeOfState))

def waveFuncMatrixElement(state1, state2):
    answer = 0
    for s1 in state1.states:
        for s2 in state2.states:
            answer += diagonalisePertubationIntegerEffect.NElectronMatrixElement(s1[1], s2[1], state1.magneticLength)*s1[0]*s2[0]
    return answer
