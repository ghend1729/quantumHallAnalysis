#Series of functions for the Slater determinent to symmetric sum polynomial bases transform
import waveFunctionClasses
import scipy
import numpy
import copy
import math

def convertSlaterToWaveFunc(slaterDet):
    return waveFunctionClasses.waveFunction([[1, slaterDet]], 1)

def S_kSingleSlat(SlaterDeterminent, k):
    answer = waveFunctionClasses.waveFunction(0,1,zeroVec=True)
    allowedMs = [m for m in SlaterDeterminent if not (m + k) in SlaterDeterminent]
    if len(allowedMs) == 0:
        answer = convertSlaterToWaveFunc(SlaterDeterminent)
    else:
        for m in allowedMs:
            x = copy.deepcopy(SlaterDeterminent)
            x.remove(m)
            newSlat = sorted([m+k] + x)
            NumFactor = ((-1)**(SlaterDeterminent.index(m) + newSlat.index(m+k)))*math.sqrt(math.factorial(m+k)/math.factorial(m))
            answer = answer + waveFunctionClasses.waveFunction([[NumFactor, newSlat]], 1)
    return answer*(2**(k/2))

def S_k(waveFunc, k):
    answer = waveFunctionClasses.waveFunction(0,1,zeroVec=True)
    for s in waveFunc.states:
        answer = answer + S_kSingleSlat(s[1], k)*s[0]
    return answer

def P_k(waveFunc, partition):
    answer = copy.deepcopy(waveFunc)
    for k in partition:
        answer = S_k(answer, k)
    return answer

def basisConversion(slaterDets, partitions, N):
    baseState = convertSlaterToWaveFunc([i for i in range(N)])
    X = numpy.array([[convertSlaterToWaveFunc(s) | P_k(baseState, p) for p in partitions] for s in slaterDets])
    print(X)
    return numpy.linalg.inv(X)