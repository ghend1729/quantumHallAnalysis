#Interaction potentials
"""
This is a library for storing functions that describe each interaction potential
as a Haldane Pseudopotential.
"""

import math
import mpmath
from usefulTools import norm2 as norm2

def exponentailRepulsion(n, magneticLength):
    return norm2(n, magneticLength)
