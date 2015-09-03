#!/usr/bin/python

# Editing: Indented with 4 spaces following the PEP8
# Classes: Imported numpy classes/module to use array function and other operator functions
# Purpose: Calculates the individual payments in Gino's ATT family plan, and compared with the 
# totalPayment: compared with indivi

import pdb

import numpy as np

def myATT_monthly_est(members, actualPayment,data,nMembers,otherFees):
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%    individual charge
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    talk = 9.99; #% do not discount, mothly charge for individual members
    if not data: #% an empty data should be the default
        data = np.array([25.0, 15.0, 0.0, 20]); #% please discount
	# pdb.set_trace()
    print '\nNOTE: Default Data: Hsiaochu($25), Ginos($15), Jinghao($0), Xia($20)\n'
    nMembers = len(data); #% Members include Hsiaochu (305-299-7934), Gino (305-299-8134), and Jinghao (305-755-7691) 

    #shared = np.around(shared_total/nMembers,decimals=2); #% please discount. (The monthly_charge of $50 in 305-2998134 account includes monthly_charge).
    if not otherFees:
        otherFees = np.array([4.48,4.48,4.48,4.48]);  #% do not discount, other_fees ($4.11) should be the default
    print 'NOTE: Default Fees: Hsiaochu($4.48), Ginos($4.48), Jinghao($4.48)\n'
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%    other variable
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    payRatio = 0.8; #% paid_ratio = discounted_cost/original_cost

    noshare = talk + data*payRatio + otherFees;
    shared = (actualPayment - np.sum(noshare,axis=0))/nMembers; # this way includes the monthly 700 min plan plus the extra fees
    individualPayments = noshare + shared;
    print shared
#    members = ['Hsiaochu','Gino','Jinghao'];
    print 'Individual Payments = talk + data*payRatio + shared + otherFees:\n'
    for i in range(nMembers): # 0 <= i < nMembers
        print '      {5} = ${0} + ${1}*{3} + ${2} + ${4} = ${6} \n' \
        .format(talk,data[i],shared,payRatio,otherFees[i],members[i],np.around(individualPayments[i],decimals=2))
    estPayment=np.sum(individualPayments,axis=0); 
    return (actualPayment,estPayment)


# Apply the function myATT_monthly_est in the main.py program
if __name__ == "__main__":
    import sys
    members = ['Hsiaochu','Gino','Jinghao','Xia'];
    actualPay = 165.26;
    nmembers = 4;
    otherFees = [4.92, 4.48, 5.15, 10.65];
    (actualPayment,estPayment) = myATT_monthly_est(members, actualPay, [], nmembers, otherFees);
    print 'Estimated Payment (${}) == Actual Payment (${}) ?\n'.format(estPayment,actualPayment) # sum over the rows (axis=0)
