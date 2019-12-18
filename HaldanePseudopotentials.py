import numpy
import mpmath
import scipy
import math
import potentials

potentialParameterFunction = potentials.coulomb

VMemory = {}

def V(m):
    if m in VMemory:
        return VMemory[m]
    else:
        answer = potentialParameterFunction(m, 1)
        VMemory[m] = answer
        return answer

def nCr(n, r):
    return math.factorial(n)//math.factorial(r)//math.factorial(n-r)

def projectToRelativeCoordinateBasis_NumpyVersion(r, s, m, M):
    sumRange = range(max([0, r - M]), min([r, m]) + 1)
    x = math.factorial(r)/math.factorial(m)
    y = math.factorial(s)/math.factorial(M)
    z = 2**(m+M)
    upFrontFactor = math.sqrt(x*y/z)
    return sum([upFrontFactor*nCr(m, a)*nCr(M, r - a)*((-1)**a) for a in sumRange])

def potential(r, s, t, u):
    if r + s == t + u:
        sumRange = range(r + s + 1)
        P = projectToRelativeCoordinateBasis_NumpyVersion
        return sum([P(r, s, m, r + s - m)*V(m)*P(t, u, m, t + u - m) for m in sumRange])
    else:
        return 0