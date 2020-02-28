from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

from sympy import *
import pprint

## Computation of the Convective Matrix
def computeA(dofs, params):
    print("\nCompute Convective Matrix \n")
    dim = params["dim"]				        # Spatial dimensions

    ## Unknown field definition
    F = DefineMatrix('F',dim+2,dim)		    # Convective Flux matrix
    Ug = dofs                               # Data interpolation to the Gauss points

    ## Other symbols definitions
    y = params["gamma"]				        # Gamma (Cp/Cv)

    ## Pgauss - Pressure definition
    pg = (y-1)*Ug[dim+1]
    for i in range(0,dim):
        pg += (y-1)*(-Ug[i+1]*Ug[i+1]/(2*Ug[0]))

    ## F - Convective Flux Matrix definition
    for j in range(0,dim):
        F[0,j] = Ug[j+1]

    for i in range (1,dim+1):
        for j in range(0,dim):
            F[i,j] = Ug[i]*Ug[j+1]/Ug[0]
            if i==j+1:
               F[i,j]+=pg

    for j in range(0,dim):
        F[dim+1,j] = (Ug[dim+1]+pg)*Ug[j+1]/Ug[0]


    ## A - Jacobian Convective Matrix definition
    A = []

    for j in range(0,dim):
        tmp = DefineMatrix('tmp',dim+2,dim+2)
        for i in range(0,dim+2):
            for k in range(0,dim+2):
      	        tmp[i,k] = diff(F[i,j], Ug[k])
      	        #print(j,'	',i,k,'=',tmp[i,k])
        A.append(tmp)
    return A

## Printing the Convective Matrix
def printA(A,params):
    dim = params["dim"]
    tmp = []
    print("The convective matrix is:\n")
    for j in range(0,dim):
    	tmp = A[j]
    	for i in range(0,dim+2):
            for k in range(0,dim+2):
                print("A[",j,",",i,",",k,"]=",tmp[i,k],"\n")

    return 0
