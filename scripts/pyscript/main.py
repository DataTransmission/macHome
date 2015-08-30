from att.plans import FamilyPlan

def info(family): 
    print family.name
    print "-------------"
    print "members: ", family.members()
    print "talk: ", family.talk
    print "total: ", family.totalPayment()
    family.printIndividualPayments()
    print "\n\n\n\n"

def main():
    family1 = FamilyPlan(
        "Gino's family",
        9.99,
        {"Gino": 25.00, "Jianhao": 15.00, "Hsiaochu": 0.00},
        {"Gino": 4.55, "Jianhao": 4.55, "Hsiaochu": 4.55},
        100.00
    )
    
    family2 = FamilyPlan(
        "Gino's family plus Doug",
        9.99,
        {"Gino": 25.00, "Jianhao": 15.00, "Hsiaochu": 0.00},
        {"Gino": 4.55, "Jianhao": 4.55, "Hsiaochu": 4.55},
        100.00
    )
    family2.data["Doug"] = 5.0
    family2.otherFees["Doug"] = 1.0
    family2.talk = 19.99
    
    #info(family1)
    #info(family2)
    family1.info()
    family2.info()

if __name__=="__main__":
    main()
