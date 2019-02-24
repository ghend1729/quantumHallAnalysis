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
#This sets the precision of all mpmath calculations. This module was used to actaully allow
#the program to acurately calculate these matrix elemnts. Essentailly these formulas have larger number
#divided by even larger number, so mpmath helps us deal with this.

def norm(m, magneticLength):
    #Calculates sqrt(<m|m>) when <z|m> = z^m*exp(-|z|^2/4*l_B). norm2 is just <m|m>
    return math.sqrt(norm2(m, magneticLength))

def allowedIJs(a, n, m):
    """
    Returns all terms in the expansion of |n>|m> in terms of relative and COM coordinates, which have
    relative angular momentum L = a.
    """
    x = ((i, j) for (i, j) in itertools.product(range(n + 1), range(m + 1)) if (i + j == a))
    return x

def IJIntigrand(a, i1, j1, i2, j2, n1, m1, n2, m2):
    """
    Retruns all combinatoric factors in the expansion of
    (Z1 - z1/2)^n1(Z1 + z1/2)^m1*(Z2 - z2/2)^n2(Z2 + z2/2)^m2
    """
    return ((-1)**(i1+i2))*nCr(n1, i1)*nCr(m1, j1)*nCr(n2, i2)*nCr(m2, j2)

def haldaneMatrixElement(n1, m1, n2, m2, magneticLength, V):
    """
    This calculates <n1|<m1|V|n2>|m2> by expanding each side in terms of relative coordinates and COM
    coordinates. V is a function V(m) which has V(m) = <m|V|m> when |m> is an unnormalised relative angular
    momentum state: <z|m> = z^m*usual exp term.
    """
    x = norm(n1, magneticLength)*norm(m1, magneticLength)*norm(n2, magneticLength)*norm(m2, magneticLength)
    if n1 + m1 == n1 + m1:
        #We need to make sure angular momentum is conserved.
        answer = 0
        for a in range(n1 + m1 + 1):
            y = V(a, magneticLength)*norm2(n1 + m1 - a, magneticLength/math.sqrt(2))
            for ((i1, j1), (i2, j2)) in itertools.product(allowedIJs(a, n1, m1), allowedIJs(a, n2, m2)):
                answer += IJIntigrand(a, i1, j1, i2, j2, n1, m1, n2, m2)*y/x/(4**a)
        return answer
    else:
        return 0
