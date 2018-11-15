#Character tables

import numpy
import scipy
import itertools
import copy
import pickle
from usefulTools import generatePartitions as partitions
import os

characterTabMemory = {}

infile = open("characterTableMemory.p", "rb")
characterTabMemory = pickle.load(infile)
infile.close()



def tabTopartition(T):
    Tconpy = copy.deepcopy(T)
    return tuple(sorted(Tconpy))

def partitionToTab(partition):
    p = copy.deepcopy(partition)
    return sorted(p, reverse=True)

def reduceTabeleau(T):
    Tconpy = copy.deepcopy(T)
    for i in range(len(Tconpy)):
        notChecked = True
        j = len(Tconpy[i]) - 1
        while notChecked:
            if j == 0:
                notChecked = False
            if Tconpy[i][j] == 0:
                Tconpy[i].pop()
            else:
                notChecked = False
            j = j - 1
    k = len(Tconpy)
    reducedRows = False
    if [] in Tconpy:
        while not reducedRows:
            k = k - 1
            if k == -1:
                reducedRows = True
            elif len(Tconpy[k]) == 0:
                Tconpy.pop()
            else:
                reducedRows = True
    return Tconpy

def checkShape(T):
    x = all([i > 0 for i in T])
    y = all([T[i+1] <= T[i] for i in range(len(T)-1)])
    return x and y

def removeHook(T, hook):
    Tconpy = copy.deepcopy(T)
    try:
        for i in range(len(T)):
            Tconpy[i] = Tconpy[i] - sum([1 for x in hook if x[0] == i])
        k = len(Tconpy)
        reducedRows = False
        if 0 in Tconpy:
            k = len(Tconpy)
            reducedRows = False
            while not reducedRows:
                k = k - 1
                if k == -1:
                    reducedRows = True
                elif Tconpy[k] == 0:
                    Tconpy.pop()
                else:
                    reducedRows = True
        return Tconpy
    except:
        print(hook)
        for i in T:
            print(i)

def genHook(startp, T, n):
    hooking = True
    hook = [startp]
    while hooking:
        if len(hook) == n:
            hooking = False
        else:
            p1 = [hook[-1][0], hook[-1][1] + 1]
            p2 = [hook[-1][0] - 1, hook[-1][1]]
            if pointInRange(p1, T) and isOnBoarder(p1, T):
                hook.append(p1)
            elif pointInRange(p2, T) and isOnBoarder(p2, T):
                hook.append(p2)
            else:
                hooking = False
    return hook

def possibleStartPoint(p, T):
    p1 = [p[0] + 1, p[1]]
    return not pointInRange(p1, T)

def pointInRange(point, T):
    if 0 <= point[0] < len(T):
        return 0 <= point[1] < T[point[0]]
    else:
        return False

def validHook(hook, T):
    return all([pointInRange(p, T) for p in hook])

def makeHooks(T, n):
    hooks = []
    startPoints = []
    for i in range(len(T)):
        for j in range(T[i]):
            startPoints.append([i,j])
    startPoints = [x for x in startPoints if possibleStartPoint(x, T)]
    hooks = [genHook(p, T, n) for p in startPoints]
    hooks = [h for h in hooks if len(h) == n]
    return hooks

def validTabsForAllHookRemoves(T, n):
    hooks = makeHooks(T, n)
    tabs = [[h, removeHook(T, h)] for h in hooks if isOnBoarderHook(h, T)]
    validTabs = [t for t in tabs if checkShape(t[1])]
    return validTabs

def isOnBoarder(p, T):
    p1 = [p[0], p[1] + 1]
    p2 = [p[0] + 1, p[1]]
    p3 = [p[0] + 1, p[1] + 1]
    x = not pointInRange(p1, T)
    y = not pointInRange(p2, T)
    z = not pointInRange(p3, T)
    return x or y or z

def isOnBoarderSpecial(p, T):
    p1 = [p[0], p[1] + 1]
    p2 = [p[0] + 1, p[1]]
    x = not pointInRange(p1, T)
    y = not pointInRange(p2, T)
    return x and y

def isOnBoarderHookEndPoint(p, T):
    p1 = [p[0], p[1] + 1]
    p2 = [p[0] + 1, p[1]]
    x = not pointInRange(p1, T)
    y = not pointInRange(p2, T)
    return x or y

def isOnBoarderStartPointSpecial(startp, secondp, T):
    if startp[1] < secondp[1]:
        p2 = [startp[0] + 1, startp[1]]
        return not pointInRange(p2, T)
    else:
        p1 = [startp[0], startp[1] + 1]
        p2 = [startp[0] + 1, startp[1]]
        x = not pointInRange(p1, T)
        y = not pointInRange(p2, T)
        return x and y

def isOnBoarderEndPointSpecial(endp, penultp, T):
    if endp[1] > penultp[1]:
        p1 = [endp[0], endp[1] + 1]
        p2 = [endp[0] + 1, endp[1]]
        x = not pointInRange(p1, T)
        y = not pointInRange(p2, T)
        return x and y
    else:
        p2 = [endp[0], endp[1] + 1]
        return not pointInRange(p2, T)

def isOnBoarderHook(h, T):
    l = len(h)
    if l == 1:
        return isOnBoarderSpecial(h[0], T)
    else:
        return isOnBoarderStartPointSpecial(h[0], h[1], T) and isOnBoarderEndPointSpecial(h[-1], h[-2], T)

def hookLen(h):
    x = max([p[0] for p in h])
    y = min([p[0] for p in h])
    return x - y

def characterTab(x, y):
    if (x, y) in characterTabMemory:
        return characterTabMemory[(x, y)]
    elif (len(x) == 0) and (len(y) == 0):
        return 1
    else:
        T = partitionToTab(x)
        y1 = y[-1]
        newy = tuple([y[i] for i in range(len(y) - 1)])
        hookedTabs = validTabsForAllHookRemoves(T, y1)
        chi = sum([((-1)**(hookLen(i[0])))*characterTab(tabTopartition(i[1]), newy) for i in hookedTabs])
        characterTabMemory[(x, y)] = chi
        return chi
"""
testHolder = 0
testPartitions = [partitions(i) for i in range(1,28)]
for i in range(27):
    for (x,y) in itertools.product(testPartitions[i], repeat=2):
        testHolder = characterTab(x,y)
    print("Done: " + str(i))

outfile = open("characterTableMemory.p", "wb")
pickle.dump(characterTabMemory, outfile)
outfile.close()

x = (1,3,4,7)
y = (1,2,3,4,5)

print(characterTab(x,y))
"""
