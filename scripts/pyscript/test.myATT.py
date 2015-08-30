#!/usr/bin/python

# Editing: Indented with 4 spaces following the PEP8
# Classes: Imported numpy classes/module to use array function and other operator functions
# Purpose: Calculates the individual payments in Gino's ATT family plan, and compared with the 
# totalPayment: compared with indivi

import pdb

import numpy as np

def myATT_monthly_est(actualPayment,data,nMembers,otherFees):
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%    700 min family plan charge 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    shared_total = 50.0; #% shared by nMembers
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%    individual charge
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    talk = 9.99; #% do not discount, mothly charge for individual members
    if not data: #% an empty data should be the default
        data = np.array([25.0,15.0,0.0]); #% please discount
	# pdb.set_trace()
    print '\nNOTE: Default Data: Hsiaochu($25), Ginos($15), Jinghao($0)\n'
    nMembers = len(data); #% Members include Hsiaochu (305-299-7934), Gino (305-299-8134), and Jinghao (305-755-7691) 

    shared = np.around(shared_total/nMembers,decimals=2); #% please discount. (The monthly_charge of $50 in 305-2998134 account includes monthly_charge).
    if not otherFees:
        otherFees = np.array([4.11,4.11,4.11]);  #% do not discount, other_fees ($4.11) should be the default
    print 'NOTE: Default Fees: Hsiaochu($4.11), Ginos($4.11), Jinghao($4.11)\n'
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%    other variable
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    payRatio = 0.8; #% paid_ratio = discounted_cost/original_cost

    individualPayments = talk + (data + shared)*payRatio + otherFees;

    members = ['Hsiaochu','Gino','Jinghao'];
    print 'Individual Payments = talk + (data + shared)*payRatio + otherFees:\n'
    for i in range(nMembers): # 0 <= i < nMembers
        print '      {5} = ${0} + (${1} + ${2})*{3} + ${4} = ${6} \n' \
        .format(talk,data[i],shared,payRatio,otherFees[i],members[i],np.around(individualPayments[i],decimals=2))
    estPayment=np.sum(individualPayments,axis=0); 
    return (actualPayment,estPayment)

def main():
    # Apply the function myATT_monthly_est in the main.py program
    # import sys
    (actualPayment,estPayment) = myATT_monthly_est(136,[],3,[]);
    print 'Estimated Payment (${}) == Actual Payment (${}) ?\n'.format(actualPayment,estPayment) # sum over the rows (axis=0)

