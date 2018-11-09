#character table alternative code perhaps more efficient
import numpy
import scipy
import itertools

from usefulTools import generatePartitions as partitions

def isCorrectShape(y):
    p1 = all(y[i] <= y[i+1] for i in range(len(y)-1))
    p2 = all([i > 0 for i in y])
    return p1 and p2

def reduce(y, x):
    reducedY = []
