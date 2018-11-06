#Character tables

import numpy
import scipy
import itertools

def tabTopartition(T):
    return tuple(sorted([len(r) for r in T]))

def reduceTabeleau(T):
    for i range(len(T)):
        notChecked = True
        j = len(T[i]) - 1
        while notChecked:
            if j == 0:
                notChecked = False
            if T[i][j] == 0:
                T[i].pop()
            else:
                notChecked = False
            j = j - 1
    return T

def checkShape(T):
    x = all(not (0 in i) for i in T)
    y = all([len(T[i+1]) <= len(T[i]) for i in range(len(T)-1)])
    return x and y

def removeHook(T, hook):
    for x in hook:
        T[x[0]][x[1]] = 0
    return reduceTabeleau(T)

def genHook(pathsAsIndicator, startp):
    hook = [startp]
    newPoint = [0,0]
    for i in pathsAsIndicator:
        if i:
            newPoint = [hook[-1][0] - 1, hook[-1][1]]
        else:
            newPoint = [hook[-1][0], hook[-1][1] + 1]
        hook.append(newPoint)

def pointInRange(point, T):
    x = point[0] < len(T)
    y = point[1] < len(T[point[0]])
    return x and y

def validHook(hook, T):
    return all([pointInRange(p, T) for p in hook])

def makeHooks(T, n):
    pathsAsIndicators = itertools.product([True, False], repeat=n)
    hooks = []
    for i in range(len(T)):
        for j in range(len(T[i])):
            tempHooks = [genHook(p, [i,j]) for p in pathsAsIndicators]
            hooks = hooks + [x for x in tempHooks if validHook(x, T)]
    return hooks

def validTabsForAllHookRemoves(T, n):
    hooks = makeHooks(T, n)
    tabs = [[h, removeHook(T, h)] for h in hooks]
    validTabs = [t for t in tabs if checkShape(t[1])]
    return validTabs
