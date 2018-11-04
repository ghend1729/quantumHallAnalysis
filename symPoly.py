#Symetric Polynomial classes
import scipy
import numpy

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

class genralSymPoly:
    def __init__(self, listOfPowerSums):
        self.powerSums = listOfPowerSums

    def __add__(self, otherObj):
        if type(otherObj) is symetricPowerSumPoly:
            listOfPartitions = [x.partition for x in self.powerSums]
            listOfPowerSums = self.powerSums
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
            for s in otherObj:
                resultPoly = resultPoly + self*s
            return resultPoly
        else:
            return genralSymPoly([s*otherObj for s in self.powerSums])

    def __pow__(self, n):
        resultPoly = self
        for i in range(n-1):
            resultPoly = resultPoly*self
        return resultPoly
