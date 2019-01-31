#Interaction potentials
"""
This is a library for storing functions that describe each interaction potential
as a Haldane Pseudopotential.
"""

import math
import scipy.special
import mpmath
from usefulTools import norm2 as norm2

def exponentailRepulsion(n, magneticLength):
    return norm2(n, magneticLength)

def v3(m, magneticLength):
    if m == 3:
        return norm2(m, magneticLength*math.sqrt(2))
    else:
        return 0

def v_k(m, magneticLength):
    if m == 9:
        return norm2(m, magneticLength*math.sqrt(2))
    else:
        return 0

def coulomb(m, magneticLength):
    return 2*math.pi*magneticLength*(2*magneticLength)**(2*m)*scipy.special.gamma(m+1/2)
