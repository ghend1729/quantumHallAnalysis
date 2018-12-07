#useful tools
import math
import mpmath
import sympy
import sympy.combinatorics
partitionMemory = {}

def signOfPermutation(x):
    p = sympy.combinatorics.permutations.Permutation(x)
    if p.is_even:
        return 1
    else:
        return -1

def generatePartitions(L):
    if L in partitionMemory:
        return partitionMemory[L]
    else:
        #source: https://stackoverflow.com/questions/10035752/elegant-python-code-for-integer-partitioning
        answer = set()
        answer.add((L, ))
        for x in range(1, L):
            for y in generatePartitions(L - x):
                answer.add(tuple(sorted((x, ) + y)))
        partitionMemory[L] = answer
    return answer

def nCr(n, r):
    f = math.factorial
    return f(n)//f(r)//f(n-r)

def norm2(n, magneticLength):
    return math.pi*mpmath.factorial(n)*(2**(n+1))*magneticLength**(2*(n+1))
