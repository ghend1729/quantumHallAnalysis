#Haldane pseudopotentials functions
"""
This module contains all functions needed to calculate 2 body matrix elements
for a particular potentail using Haldane pseudopotentials to do so.
"""

import numpy
import math
import mpmath
import itertools
from usefulTools import norm2 as norm2
from usefulTools import nCr as nCr
mpmath.mp.dps = 25

def norm(m, magneticLength):
    return math.sqrt(norm2(m, magneticLength))

def allowedIJs(a, n, m):
    x = ((i, j) for (i, j) in itertools.product(range(n + 1), range(m + 1)) if (i + j == a))
    return x

def IJIntigrand(a, i1, j1, i2, j2, n1, m1, n2, m2):
    return ((-1)**(i1+i2))*nCr(n1, i1)*nCr(m1, j1)*nCr(n2, i2)*nCr(m2, j2)

def haldaneMatrixElement(n1, m1, n2, m2, magneticLength, V):
    x = norm(n1, magneticLength)*norm(m1, magneticLength)*norm(n2, magneticLength)*norm(m2, magneticLength)
    if n1 + m1 == n1 + m1:
        answer = 0
        for a in range(n1 + m1 + 1):
            y = V(a, magneticLength)*norm2(n1 + m1 - a, magneticLength/math.sqrt(2))
            for ((i1, j1), (i2, j2)) in itertools.product(allowedIJs(a, n1, m1), allowedIJs(a, n2, m2)):
                answer += IJIntigrand(a, i1, j1, i2, j2, n1, m1, n2, m2)*y/x/(4**a)
        return answer
    else:
        return 0
