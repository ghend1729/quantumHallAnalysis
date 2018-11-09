#Symetric Polynomial classes
import scipy
import numpy
import characterTables
from usefulTools import generatePartitions

class symetricPowerSumPoly:
    def __init__(self, partition, coeficient):
        self.partition = tuple(sorted(partition))
        self.coeficient = coeficient

    def __add__(self, otherPoly):
        if type(otherPoly) is symetricPowerSumPoly:
            if self.partition == otherPoly.partition:
                return symetricPowerSumPoly(self.partition, self.coeficient + otherPoly.coeficient)
            else:
                return genralSymPoly([self, otherPoly])
        else:
            return otherPoly + self

    def __mul__(self, otherObj):
        if type(otherObj) is symetricPowerSumPoly:
            return symetricPowerSumPoly(self.partition + otherObj.partition, self.coeficient*otherObj.coeficient)
        elif type(otherObj) is genralSymPoly:
            return otherObj*self
        else:
            return symetricPowerSumPoly(self.partition, self.coeficient*otherObj)

    def convertToShur(self):
        partitions = generatePartitions(sum(self.partition))
        resultPoly = genralSymPolySchur([schurPoly(p, self.coeficient*characterTables.characterTab(p, self.partition)) for p in partitions])
        return resultPoly

    def __str__(self):
        return str(self.coeficient) + "*" + str(self.partition)

class schurPoly:
    def __init__(self, partition, coeficient):
        self.partition = tuple(sorted(partition))
        self.coeficient = coeficient

    def __add__(self, otherPoly):
        if type(otherPoly) is schurPoly:
            if self.partition == otherPoly.partition:
                return schurPoly(self.partition, self.coeficient + otherPoly.coeficient)
            else:
                return genralSymPolySchur([self, otherPoly])
        else:
            return otherPoly + self

    def __mul__(self, otherObj):
        return schurPoly(self.partition, self.coeficient*otherObj)

    def __str__(self):
        return str(self.coeficient) + "*" + str(self.partition)

class genralSymPoly:
    def __init__(self, listOfPowerSums):
        self.powerSums = [s for s in listOfPowerSums if not (s.coeficient == 0)]

    def __add__(self, otherObj):
        if type(otherObj) is symetricPowerSumPoly:
            listOfPartitions = [x.partition for x in self.powerSums]
            listOfPowerSums = [x for x in self.powerSums]
            if otherObj.partition in listOfPartitions:
                for i in range(len(listOfPowerSums)):
                    if listOfPowerSums[i].partition == otherObj.partition:
                        listOfPowerSums[i] = listOfPowerSums[i] + otherObj
                return genralSymPoly(listOfPowerSums)
            else:
                return genralSymPoly(listOfPowerSums + [otherObj])
        else:
            resultPoly = self
            for s in otherObj.powerSums:
                resultPoly = resultPoly + s
            return resultPoly

    def __mul__(self, otherObj):
        if type(otherObj) is genralSymPoly:
            resultPoly = genralSymPoly([])
            for s in otherObj.powerSums:
                resultPoly = resultPoly + self*s
            return resultPoly
        else:
            return genralSymPoly([s*otherObj for s in self.powerSums])

    def __pow__(self, n):
        resultPoly = self
        for i in range(n-1):
            resultPoly = resultPoly*self
        return resultPoly

    def convertToShurBasis(self):
        resultPoly = genralSymPolySchur([])
        for s in self.powerSums:
            resultPoly = resultPoly + s.convertToShur()
        return resultPoly

    def __str__(self):
        resultString = ""
        for i in range(len(self.powerSums)-1):
            resultString += str(self.powerSums[i]) + " + "
        resultString += str(self.powerSums[-1])
        return resultString

class genralSymPolySchur:
    def __init__(self, listOfShurs):
        self.schurs = [s for s in listOfShurs if not (s.coeficient == 0)]

    def __add__(self, otherObj):
        if type(otherObj) is schurPoly:
            listOfPartitions = [x.partition for x in self.schurs]
            listOfSchurs = [x for x in self.schurs]
            if otherObj.partition in listOfPartitions:
                for i in range(len(listOfSchurs)):
                    if listOfSchurs[i].partition == otherObj.partition:
                        listOfSchurs[i] = listOfSchurs[i] + otherObj
                return genralSymPolySchur(listOfSchurs)
            else:
                return genralSymPolySchur(listOfSchurs + [otherObj])
        else:
            resultPoly = self
            for s in otherObj.schurs:
                resultPoly = resultPoly + s
            return resultPoly

    def __mul__(self, otherObj):
        return genralSymPolySchur([s*otherObj for s in self.schurs])

    def __str__(self):
        resultString = ""
        for i in range(len(self.schurs)-1):
            resultString += str(self.schurs[i]) + " + "
        resultString += str(self.schurs[-1])
        return resultString

z = symetricPowerSumPoly((1,2), 1)
print(z.convertToShur())
x = symetricPowerSumPoly((1,1,1), 1)
y = z + x
print(y.convertToShurBasis())
print((y**2).convertToShurBasis())
