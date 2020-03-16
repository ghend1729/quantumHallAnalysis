import numpy
import mpmath
import scipy
import math
import mpmath
import potentials

useMpmath = True

mpmath.mp.dps = 30

potentialParameterFunction = potentials.longRange

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

def nCrMP(n, r):
    f = mpmath.factorial
    return f(n)/f(r)/f(n-r)

def projectToRelativeCoordinateBasis_NumpyVersion(r, s, m, M):
    sumRange = range(max([0, r - M]), min([r, m]) + 1)
    x = math.factorial(r)/math.factorial(m)
    y = math.factorial(s)/math.factorial(M)
    z = 2**(m+M)
    upFrontFactor = math.sqrt(x*y/z)
    return sum([upFrontFactor*nCr(m, a)*nCr(M, r - a)*((-1)**(m-a)) for a in sumRange])

def projectToRelativeCoordinateBasis_mpmathVersion(r, s, m, M):
    sumRange = range(max([0, r - M]), min([r, m]) + 1)
    x = mpmath.factorial(r)/mpmath.factorial(m)
    y = mpmath.factorial(s)/mpmath.factorial(M)
    z = 2**(m+M)
    upFrontFactor = mpmath.sqrt(x*y/z)
    return sum([upFrontFactor*nCrMP(m, a)*nCrMP(M, r - a)*((-1)**a) for a in sumRange])

def projectToRelativeCoordinateBasis(r, s, m, M):
    if useMpmath:
        return projectToRelativeCoordinateBasis_mpmathVersion(r, s, m, M)
    else:
        return projectToRelativeCoordinateBasis_NumpyVersion(r, s, m, M)

def potential(magneticLength, r, s, t, u):
    if r + s == t + u:
        sumRange = range(1, r + s + 1, 2)
        P = projectToRelativeCoordinateBasis
        return sum([P(r, s, m, r + s - m)*V(m)*P(t, u, m, t + u - m) for m in sumRange])
    else:
        return 0