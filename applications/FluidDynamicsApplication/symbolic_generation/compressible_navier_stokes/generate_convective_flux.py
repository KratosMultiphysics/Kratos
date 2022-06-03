from sympy import *
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

def ComputeEulerJacobianMatrix(dofs, params):
    """This function calculates the Euler Jacobian matrix for convection"""

    print("\nCompute Euler Jacobian matrix \n")

    ## Auxiliary variables
    dim = params["dim"]
    gamma = params["gamma"]
    rho = dofs[0]
    mom = []
    vel = []
    mom_prod = 0.0
    for i in range(dim):
        mom.append(dofs[i + 1])
        vel.append(dofs[i + 1] / rho)
        mom_prod += dofs[i + 1]**2
    e_tot = dofs[dim + 1]
    p = (gamma - 1) * (e_tot - 0.5 * mom_prod / rho)

    ## Define and fill the convective flux matrix
    E = DefineMatrix('E', dim + 2, dim)
    for j in range(dim):
        E[0,j] = mom[j]
        for d in range(dim):
            E[1 + d,j] = mom[j]*vel[d]
            if j == d:
                E[1 + d,j] += p
        E[dim + 1,j] = vel[j] * (e_tot + p)

    ## Define and derive the Euler Jacobian matrix
    A = []
    for j in range(dim):
        A_j = DefineMatrix('A_j', dim + 2, dim + 2)
        for m in range(dim + 2):
            for n in range(dim + 2):
                A_j[m,n] = diff(E[m,j], dofs[n])
        A.append(A_j)

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
