"""
These are functions for elementary tests in criteria.
They accept either two or three argument, depending on the operation.
The output argument is a boolean.
"""
def isLowerThanOrEqualTo(a,b):
    return a<=b

def isStrictlyLowerThan(a,b):
    return a<b

def isGreaterThanOrEqualTo(a,b):
    return a>=b

def isStrictlyGreaterThan(a,b):
    return a>b

def isDistanceLowerThanOrEqualTo(a,b,tolerance):
    return isLowerThanOrEqualTo(abs(a-b),tolerance)

def isDistanceStrictlyLowerThan(a,b,tolerance):
    return isStrictlyLowerThan(abs(a-b),tolerance)

def isDifferenceLowerThanOrEqualTo(a,b,tolerance):
    return isLowerThanOrEqualTo(a-b,tolerance)

def isDifferenceStrictlyLowerThan(a,b,tolerance):
    return isStrictylLowerThan(a-b,tolerance)

def isEqualTo(a,b):
    return a==b
