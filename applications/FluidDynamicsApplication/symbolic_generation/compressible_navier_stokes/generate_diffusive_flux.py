from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

from sympy import *
import pprint

## Computation of the Diffusive Matrix
def computeG(dofs,params,Hg,Gg):
    print("\nCompute Diffusive Matrix\n")
    dim = params["dim"]				                    # spatial dimensions

    ## Unknown fields definition
    H = Hg.copy()                               		# Gradient of U
    G = Gg.copy()                               		# Diffusive Flux matrix
    tau_stress = DefineMatrix('tau_stress',dim,dim)		# Shear stress tensor for Newtonian fluid
    q = DefineVector('q',dim)			                # Heat flux vector

    ## Other simbols definition
    c_v = params["c_v"]				                    # Specific Heat at Constant volume
    gamma = params["gamma"]				                # Gamma (Cp/Cv)
    mu  = params["mu"]         			                # Dynamic viscosity
    l = params["lambda"]			                    # Thermal Conductivity of the fluid

    ## Data interpolation to the Gauss points
    Ug = dofs

    ## Pgauss - Pressure definition
    pg = Ug[dim+1]
    for i in range(0,dim):
        pg += (-Ug[i+1]*Ug[i+1]/(2*Ug[0]))
    pg *= (gamma-1)

    ## Tau - Shear stress tensor definition
    for i in range(0,dim):
        for j in range(i,dim):
            if i!=j:
               tau_stress[i,j] = (mu/Ug[0])*(H[i+1,j]+H[j+1,i])-(mu/Ug[0]**2)*(Ug[i+1]*H[0,j]+Ug[j+1]*H[0,i])
            if i==j:
               tau_stress[i,j]= (2*mu/Ug[0])*H[i+1,i]-(2*mu/Ug[0]**2)*Ug[i+1]*H[0,i]
               for k in range(0,dim):
                   tau_stress[i,j]+= -(2*mu/(3*Ug[0]))*H[k+1,k]+(2*mu/(3*Ug[0]**2))*Ug[k+1]*H[0,k]

    for i in range(1,dim):
        for j in range(0,dim-1):
            if j!=i:
               tau_stress[i,j] = tau_stress[j,i]

    ## q - Heat flux vector definition
    for i in range(0,dim):
        q[i] = l*Ug[dim+1]/(Ug[0]**2*c_v)*H[0,i]-(l*H[dim+1,i])/(Ug[0]*c_v)
        for j in range(0,dim):
            q[i] += -l*Ug[j+1]**2/(c_v*Ug[0]**3)*H[0,i]+l/(Ug[0]**2*c_v)*Ug[j+1]*H[j+1,i]
    #NB!!!There is an error in the definition of q[i] in the research proposal.
    #The second term of the equation has an opposite sign!!!NB#
    '''
    G [(dim+2)*(dim)]

    0                                   0
    -tau00                              -tau01
    -tau01                              -tau11
    -mu/rho*tau00-mv/rho*tau01+q0       -mu/rho*tau01-mv/rho*tau11+q1
    '''
    ## G - Diffusive Matrix definition
    for j in range(0,dim):
        G[0,j]= 0 			                            #Mass equation related

    for i in range(1,dim+1):
        for j in range(0,dim):
            G[i,j]=-tau_stress[i-1,j]		            #Moment equation related

    for j in range(0,dim):                              #Energy equation related
        G[dim+1,j] = q[j]
        for k in range(0,dim):
            G[dim+1,j] += -Ug[k+1]*tau_stress[k,j]/Ug[0]

    return G


## Computation of the Diffusive Matrix with Shock Capturing
def computeGsc(dofs,params,Hg,Gg,v_sc,k_sc):
    print("\nCompute Diffusive Matrix with Shock Capturing\n")
    dim = params["dim"]				                    # spatial dimensions

    ## Unknown fields definition
    H = Hg.copy()                               		# Gradient of U
    Gsc = Gg.copy()                               		# Diffusive Flux matrix
    tau_stress = DefineMatrix('tau_stress',dim,dim)		# Shear stress tensor for Newtonian fluid
    q = DefineVector('q',dim)			                # Heat flux vector

    ## Other simbols definition
    c_v = params["c_v"]				                    # Specific Heat at Constant volume
    gamma = params["gamma"]				                # Gamma (Cp/Cv)
    mu  = params["mu"]         			                # Dynamic viscosity
    l = params["lambda"]			                    # Thermal Conductivity

    ## Data interpolation to the Gauss points
    Ug = dofs

    ## Pgauss - Pressure definition
    pg = Ug[dim+1]
    for i in range(0,dim):
        pg += (-Ug[i+1]*Ug[i+1]/(2*Ug[0]))
    pg *= (gamma-1)

    ## tau - Shear stress tensor definition
    for i in range(0,dim):
        for j in range(i,dim):
            if i!=j:
               tau_stress[i,j] = (mu/Ug[0])*(H[i+1,j]+H[j+1,i])-(mu/Ug[0]**2)*(Ug[i+1]*H[0,j]+Ug[j+1]*H[0,i])
            if i==j:
               tau_stress[i,j]= (2*mu/Ug[0])*H[i+1,i]-(2*mu/Ug[0]**2)*Ug[i+1]*H[0,i]
               for k in range(0,dim):
                   tau_stress[i,j]+= -(2*mu/(3*Ug[0]))*H[k+1,k]+(2*mu/(3*Ug[0]**2))*Ug[k+1]*H[0,k]

    for i in range(1,dim):
        for j in range(0,dim-1):
            if j!=i:
               tau_stress[i,j] = tau_stress[j,i]

    ## q - Heat flux vector definition
    for i in range(0,dim):
        q[i] = l*Ug[dim+1]/(Ug[0]**2*c_v)*H[0,i]-(l*H[dim+1,i])/(Ug[0]*c_v)
        for j in range(0,dim):
            q[i] += -l*Ug[j+1]**2/(c_v*Ug[0]**3)*H[0,i]+l/(Ug[0]**2*c_v)*Ug[j+1]*H[j+1,i]
    #NB!!!There is an error in the definition of q[i] in the research proposal.
    #The second term of the equation has an opposite sign!!!NB#
    '''
    G [(dim+2)*(dim)]

    0                                   0
    -tau00                              -tau01
    -tau01                              -tau11
    -mu/rho*tau00-mv/rho*tau01+q0       -mu/rho*tau01-mv/rho*tau11+q1
    '''
    tau_sc = (1+(Ug[0]*v_sc)/mu)*tau_stress     # Stress tensor with shock capturing viscosity
    q_sc = (1+(Ug[0]*c_v*k_sc)/l)*q             # Heat flux with shock capturing conductivity

    ## Gsc - Diffusive Matrix definition
    for j in range(0,dim):
        Gsc[0,j]= 0 			                #Mass equation related

    for i in range(1,dim+1):
        for j in range(0,dim):
            Gsc[i,j]=-tau_sc[i-1,j]		        #Moment equation related

    for j in range(0,dim):                      #Energy equation related
        Gsc[dim+1,j] = q_sc[j]
        for k in range(0,dim):
            Gsc[dim+1,j] += -Ug[k+1]*tau_sc[k,j]/Ug[0]

    return Gsc

## Printing the Diffusive Matrix G
def printK(G,params):
    dim = params["dim"]
    print("The diffusive matrix is:\n")
    for ll in range(dim+2):
        for mm in range(dim):
            print("G[",ll,",",mm,"]=",G[ll,mm],"\n")

    return 0
