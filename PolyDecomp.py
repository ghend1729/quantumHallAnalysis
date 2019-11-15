import itertools
import usefulTools
import numpy
import math

class polyNom:
    def __init__(self, polyList):
        self.polyL = polyList

    def __mul__(self, x):
        newList = [[a[0]*b[0], [a[1][i] + b[1][i] for i in range(len(a[1]))]] for (a, b) in itertools.product(self.polyL, x.polyL)]
        return polyNom(newList)

    def dot(self, x):
        return sum([a[0]*b[0]*kDelta(a[1], b[1]) for (a, b) in itertools.product(self.polyL, x.polyL)])



def kDelta(a, b):
    if all([a[i] == b[i] for i in range(len(a))]):
        return 1
    else:
        return 0

def makeSymSum(a, N):
    polyList = [[1, list(x)] for x in set(itertools.permutations([a] + (N-1)*[0]))]
    return polyNom(polyList)

def genSlaterDet(a, N):
    bareList = [i for i in range(N)]
    polyList = [[usefulTools.signOfPermutation(x),[a[x[i]] for i in range(N)]] for x in itertools.permutations(bareList)]
    return polyNom(polyList)

def antiSymPolySum(partition, N, magneticLength):
    bareList = [i for i in range(N)]
    finalPoly = genSlaterDet(bareList)
    for p in partition:
        finalPoly = makeSymSum(p, N)*finalPoly
    finalPoly.polyL = [[a[0]*numpy.product([math.sqrt(usefulTools.norm2(x, magneticLength)) for x in a[1]]), a[1]] for a in finalPoly.polyL]
    return finalPoly

