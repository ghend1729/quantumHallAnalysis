#Coulomb matrix element functions
#These functions are for calculating relavent matrix elements for the Coulomb interaction between TWO electrons in the lowest
#Landau level.
#
#These formulea are taken directly from 'Composite Fermions', by J. K. Jain.
#The names of the functions are the symbols given in the book for consistencyself.

import math
import scipy
import scipy.special
import mpmath

mpmath.mp.dps = 50

def AIntegrand(r,s,t,i):
    f = mpmath.gamma
    g = math.factorial
    x = mpmath.binomial(s,i)*(f(i+1/2)/g(r+i))*(f(1/2+r+i)/f(3/2+r+t+i))
    return x

def BIntegrand(r,s,t,i):
    return AIntegrand(r,s,t,i)*(1/2+r+2*i)

def A(r,s,t):
    return sum([AIntegrand(r,s,t,i) for i in range(s+1)])

def B(r,s,t):
    return sum([BIntegrand(r,s,t,i) for i in range(s+1)])

def matrixElement(magneticLength, m1Prime, m2Prime, m1, m2):
    if (m1Prime + m2Prime == m1 + m2):
        if m2 >= m2Prime:
            t = m2Prime
            s = m1
            r = m2 - t
            f = math.factorial
            g = mpmath.gamma

            x = mpmath.sqrt((f(s+r)/f(s))*(f(t+r)/f(t)))
            y = g(r+s+t+3/2)/(mpmath.pi*mpmath.power(2, r+s+t+2))
            z = A(r,s,t)*B(r,t,s) + B(r,s,t)*A(r,t,s)

            return x*y*z/magneticLength
        else:
            return matrixElement(magneticLength,m1,m2,m1Prime,m2Prime)
    else:
        return 0
