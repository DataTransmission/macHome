#!/usr/bin/python

#import pdb

#import numpy as np

class FamilyPlan:
    def __init__(self, name, talk, data, otherFees, sharedTotal):
        self.name = name
        self.talk = talk
        self.data = data
        self.otherFees = otherFees
        self.sharedTotal = sharedTotal

    def members(self):
	    return self.data.keys()

    def individualPayment(self, member):
	    return self.talk + (self.data[member] + self.sharedTotal / len(self.members()))*0.8 + self.otherFees[member]
	
    def totalPayment(self):
        total = 0.0
        for member in self.members():
            total += self.individualPayment(member)
	    return total

    def printIndividualPayments(self):
	    for member in self.members():
	        print "Member: %s Payment: %f" % (member, self.individualPayment(member))

    def info(self): 
        print self.name
        print "-------------"
        print "members: ", self.members()
        print "talk: ", self.talk
        print "total: ", self.totalPayment()
        self.printIndividualPayments()
        print "\n\n\n\n"
