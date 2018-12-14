#CFT matrix elements
import math
import itertools

def ladderOperatorMatrixElements(ups, downs, gamma, state1, state2):
    answer = 1
    f = math.factorial
    ladders = tuple([-i for i in ups]) + downs
    for i in set(ups):
        occupancyNumber = sum([1 for j in state1 if j == i])
        numberToTakeAway = sum([1 for j in ups if j == i])
        answer = answer*math.sqrt(f(occupancyNumber)//f(occupancyNumber - numberToTakeAway))
    answer = answer*sum([operatorPrefactor(x, gamma) for x in set(itertools.permutations(ladders))])
    return answer

def operatorPrefactor(ladders, gamma):
    answer = 1
    for i in range(len(ladders)):
        answer = answer*math.sqrt(abs(ladders[i]))*kappa(ladders[i], gamma[i])
    return answer

def kappa(n, gamma):
    answer = 1
    for i in range(gamma - 1):
        answer = answer*(n + i + 1)
    return answer

def TOperatorMatrixElement(gamma, state1, state2):
    f = makeLadderOperator
    if len(gamma) > len(state1) + len(state2):
        return 0
    else:
        f = makeLadderOperator
        answer = 0
        for i in (j for j in range(1, len(gamma)) if len(state1) - j == len(state2) - len(gamma) + j and len(state1) - j >= 0 and len(state2) - len(gamma) + j >= 0):
            for x in itertools.combinations(range(len(state1)), i):
                answer += sum([ladderOperatorMatrixElements(f(state1, x), f(state2, y), gamma, state1, state2) for y in itertools.combinations(range(len(state2)), len(gamma) - i) if checkNonZeroTerm(y, x, state2, state1)])
        return answer*((-1)**len(gamma))

def checkNonZeroTerm(downsIndex, upsIndex, state2, state1):
    return tuple([state1[i] for i in range(len(state1)) if not i in downsIndex]) == tuple([state2[i] for i in range(len(state2)) if not i in upsIndex])

def makeLadderOperator(state, ladderIndex):
    return tuple([state[i] for i in ladderIndex])
