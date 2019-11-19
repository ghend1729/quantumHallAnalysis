#Interaction potentials
"""
This is a library for storing functions that describe each interaction potential
as Haldane pseudo-potentials.
"""

import math
import scipy.special
import mpmath
from usefulTools import norm2 as norm2
import matplotlib
from matplotlib import pyplot as pyplot

V_m = [math.log(i+1) for i in range(10)]
mMax = len(V_m)

def exponentailRepulsion(n, magneticLength):
    return norm2(n, magneticLength)

def v3(m, magneticLength):
    if m == 3:
        return norm2(m, magneticLength*math.sqrt(2))
    else:
        return 0

def v_k(m, magneticLength):
    if m < mMax:
        return V_m[m]*norm2(m, magneticLength*0.1)
    else:
        return 0
    

def coulomb(m, magneticLength):
    return (2*math.pi)*magneticLength*((2*magneticLength)**(2*m))*scipy.special.gamma(m+1/2)

def longRange(m, magneticLength):
    alpha = -2
    return -0.5*2*math.pi*magneticLength*(2*magneticLength)**(2*m+1-alpha)*scipy.special.gamma(m+1-(alpha/2))
"""
y = [1-longRange(i, 1)/norm2(i, math.sqrt(2)) for i in range(30)]
x = [i for i in range(30)]
pyplot.plot(x,y)
pyplot.show()
"""