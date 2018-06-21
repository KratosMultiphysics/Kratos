# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
from KratosMultiphysics import *
import random
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.EigenSolversApplication import *
from custom_timer import Timer

size = 10000000
size_c = 2
M = Matrix(size_c, size)

timer = Timer()
timer.StartTimer()
timer.StartNewLap()

#for i in range(size_c):
#    for j in range(size):
#        M[i,j] = random.uniform(-10.0, 100.0)
for i in range(size):
#    for j in range(size):
    M[0,i] = random.uniform(-10.0, 100.0)
    M[1,i] = random.uniform(1000,2000)
print("\n> Time needed for assign loop = ", timer.GetLapTime(), "s")
U = Matrix()
V = Matrix()
S = Matrix()

svd_utility = SingularValueDecomposition()


timer.StartTimer()
timer.StartNewLap()
print("\n AlgorithmSVD: SVD start!!")
svd_utility.Solve(M,U,V,S)
print("\n AlgorithmSVD: SVD finished!!")
print("\n> Time needed for SVD = ", timer.GetLapTime(), "s")
#print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")
#print(U)
#print(V.size_1)
#print(V.size_2)