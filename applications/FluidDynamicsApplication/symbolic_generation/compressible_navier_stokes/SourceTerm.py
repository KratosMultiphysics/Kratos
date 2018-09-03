from KratosMultiphysics import *

from sympy import *
from sympy_fe_utilities import *
import pprint

## Computation of the Source Matrix
def computeS(force,source,params):
    print("\nCompute Source Matrix \n")
    dim = params["dim"]				        # Spatial dimensions
    
    ## Unknown field definition	
    f = force.copy()					    # Body force vector
    r = source			 		            # Heat Source/Sink Term  #Symbol('r', positive = True) 
    '''
    S - Reactive Matrix definition     
     0  0  0   0
     fx 0  0   0
     fy 0  0   0
     r  fx fy  0
    '''
    S = zeros(dim+2,dim+2)                  # Reactive matrix (Source terms)

    
    for i in range(1,dim+1):
        S[i,0]   = f[i-1]
        S[dim+1,i] = f[i-1]        
    S[dim+1,0] = r
        
    return S

## Printing of the Source Matrix   
def printS(S,params):
    dim = params["dim"]				        # Spatial dimensions
    print("The source term matrix is:\n")
    for i in range (0,dim+2):
        for j in range (0,dim+2):
            print("S[",i,",",j,"]=",S[i,j],"\n")

    return 0
