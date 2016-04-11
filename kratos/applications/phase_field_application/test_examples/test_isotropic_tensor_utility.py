import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.PhaseFieldApplication import *
kernel = Kernel()   #defining kernel

tester = IsotropicTensorUtilityTester()
# tester.Test1()
# print("End Test 1----------------------------------")
# tester.Test1_1()
# print("End Test 1-1----------------------------------")
# tester.Test1_2()
# print("End Test 1-2----------------------------------")
tester.Test2()
print("End Test 2-----------------------------------------------------------------------")
tester.Test2_1()
print("End Test 2-1---------------------------------------------------------------------")
tester.Test2_1()
print("End Test 2-2---------------------------------------------------------------------")

