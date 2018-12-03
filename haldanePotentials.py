#Haldane pseudopotentials functions
"""
This module contains all functions needed to calculate 2 body matrix elements
for a particular potentail using Haldane pseudopotentials to do so.
"""

import numpy
import itertools

def allowedBs(a, n, m):
    x = ((i, j) for (i, j) in itertools.product(range(n + 1), range(m + 1)) if (i + j == a))
    return (i - j for (i, j) in x)
    
