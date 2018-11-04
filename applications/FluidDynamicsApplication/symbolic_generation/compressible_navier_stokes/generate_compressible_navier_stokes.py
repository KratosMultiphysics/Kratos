from KratosMultiphysics import *

from sympy import *
from sympy_fe_utilities import *
import pprint

from params_dict import params
import ConvectiveFlux
import DiffusiveFlux
import SourceTerm
import StabilizationMatrix


dim = params["dim"]                         # Define Dimension in params.py  
BlockSize = dim+2 					        # Dimension of the vector of Unknowns
do_simplifications = False
mode = "c"                                  # Output mode to a c++ file


if(dim == 2):
    nnodes = 3
elif(dim == 3):
    nnodes = 4

impose_partion_of_unity = False
N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

# Unknown fields definition (Used later for the gauss point interpolation)
U = DefineMatrix('U',nnodes,BlockSize)	     # Vector of Unknowns ( Density,Velocity[dim],Total Energy )
Un = DefineMatrix('Un',nnodes,BlockSize)     # Vector of Unknowns one step back
Unn = DefineMatrix('Unn',nnodes,BlockSize)   # Vector of Unknowns two steps back
r = DefineVector('r',nnodes)                 # Sink term    #COMMENT for manufactured solution

# Test functions defintiion
w = DefineMatrix('w',nnodes,BlockSize)	    # Variables field test

# External terms definition
f_ext = DefineMatrix('f_ext',nnodes,dim)    # Forcing term #COMMENT for manufactured solution

# Definition of other symbols
bdf0 = Symbol('bdf0')                       # Backward differantiation coefficients
bdf1 = Symbol('bdf1')
bdf2 = Symbol('bdf2')
v_sc = Symbol('v_sc')                       # Shock capturing Viscosity
k_sc = Symbol('k_sc')                       # Shock capturing Conductivity

### Construction of the variational equation
Ug = DefineVector('Ug',BlockSize)			# Dofs vector
H = DefineMatrix('H',BlockSize,dim)			# Gradient of U
f = DefineVector('f',dim)			        # Body force vector
rg = Symbol('rg', positive = True)		    # Source/Sink term
V = DefineVector('V',BlockSize)			    # Test function
Q = DefineMatrix('Q',BlockSize,dim)			# Gradient of V
acc = DefineVector('acc',BlockSize)         # Derivative of Dofs/Time
#G = DefineMatrix('G',BlockSize,dim)		# Diffusive Flux matrix
Gsc = DefineMatrix('G',BlockSize,dim)       # Diffusive Flux matrix with Shock Capturing


## Matrix Computation

S = SourceTerm.computeS(f,rg,params)
#SourceTerm.printS(S,params)
A = ConvectiveFlux.computeA(Ug,params)
#ConvectiveFlux.printA(A,params)
#G = DiffusiveFlux.computeG(Ug,params,H,G)
Gsc = DiffusiveFlux.computeGsc(Ug,params,H,Gsc,v_sc,k_sc)   
#DiffusiveFlux.printK(Gsc,params)
Tau = StabilizationMatrix.computeTau(params)
#StabilizationMatrix.printTau(Tau,params)

## Nonlinear operator definition   
l1 = Matrix(zeros(dim+2,1))		            # Convective Matrix*Gradient of U
A_small = []
for j in range(0,dim):
    A_small = A[j]
    for ll in range(BlockSize):
        for mm in range(BlockSize):
            l1[ll] += A_small[ll,mm]*H[mm,j]

l3 = S*Ug				                    # Source term
print("\nCompute Non-linear operator\n")
L = l1-l3                                   # Nonlinear operator

## Residual definition     
res = -acc - L		

## Nonlinear adjoint operator definition  
m1 = Matrix(zeros(dim+2,1))		            # Convective term
psi = Matrix(zeros(dim+2,dim))

for j in range(0,dim):
    A_T = A[j].transpose()
    for l in range(0,dim+2):
        for m in range(0,dim+2):
            psi[l,j] += A_T[l,m]*Q[m,j]                 
            for n in range(0,dim+2):
                psi[l,j] +=diff(A_T[l,m],Ug[n])*H[n,j]*V[m]   

for s in range(0,dim+2):
    for j in range(0,dim):
        m1[s] += psi[s,j]

m3 = S.transpose()*V			            # Source term

L_adj = -m1-m3                              # Nonlinear adjoint operator

## Istotropic Residual Based Shock Capturing
res_m = Matrix(zeros(dim,1))                # Momentum residual
for i in range(0,dim):
    res_m[i,0] = res[i+1,0]

res_e = Matrix(zeros(1,1))                  # Energy residual
res_e[0,0] = res[dim+1]


## Variational Formulation - Final equation
n1 = V.transpose()*acc		                # Mass term - FE scale
    
temp = zeros(dim+2,1)
A_smalll = []
for i in range(0,dim):
    A_smalll = A[i]
    for ll in range(BlockSize):
        for mm in range(BlockSize):
            temp[ll] += A_smalll[ll,mm]*H[mm,i]

n2 = V.transpose()*temp			            # Convective term - FE scale

n3 = Matrix(zeros(1,1))                     # Diffusive term - FE scale

for j in range(0,dim):
    for k in range(BlockSize):
        n3[0,0] += Q[k,j]*(-Gsc[k,j])       # G with shock capturing - FE scale

n4 = -V.transpose()*(S*Ug)		            # Source term - FE scale

n5 = L_adj.transpose()*(Tau*res)	        # VMS_adjoint - Subscales 

print("\nCompute Variational Formulation\n")
rv = n1+n2+n3+n4+n5 			            # VARIATIONAL FORMULATION - FINAL EQUATION

### Substitution of the discretized values at the gauss points
print("\nSubstitution of the discretized values at the gauss points\n")

## Data interpolation at the gauss points
U_gauss = U.transpose()*N
w_gauss = w.transpose()*N
f_gauss = f_ext.transpose()*N                     #COMMENT for manufactured solution
acc_gauss = (bdf0*U+bdf1*Un+bdf2*Unn).transpose()*N
r_gauss = (r.transpose()*N)[0]                    #COMMENT for manufactured solution  
#r_gauss = Symbol('r_gauss', positive = True)     #USED fro manufactured solution

## Gradients computation
grad_U = DfjDxi(DN,U).transpose()
grad_w = DfjDxi(DN,w).transpose()

print("\nSubstitution in the variational formulation\n")
SubstituteMatrixValue(rv, Ug, U_gauss)
SubstituteMatrixValue(rv, acc, acc_gauss)
SubstituteMatrixValue(rv, H, grad_U)
SubstituteMatrixValue(rv, V, w_gauss)
SubstituteMatrixValue(rv, Q, grad_w)
SubstituteMatrixValue(rv, f, f_gauss)           #COMMENT for manufactured solution
SubstituteScalarValue(rv, rg, r_gauss)          #COMMENT for manufactured solution

print("\nSubstitution in the residual of momentum\n")
SubstituteMatrixValue(res_m, Ug, U_gauss)
SubstituteMatrixValue(res_m, acc, acc_gauss)
SubstituteMatrixValue(res_m, H, grad_U)
SubstituteMatrixValue(res_m, f, f_gauss)       #COMMENT for manufactured solution
SubstituteScalarValue(res_m, rg, r_gauss)      #COMMENT for manufactured solution

print("\nSubstitution in the residual of total energy\n")
SubstituteMatrixValue(res_e, Ug, U_gauss)
SubstituteMatrixValue(res_e, acc, acc_gauss)
SubstituteMatrixValue(res_e, H, grad_U)
SubstituteMatrixValue(res_e, f, f_gauss)       #COMMENT for manufactured solution
SubstituteScalarValue(res_e, rg, r_gauss)      #COMMENT for manufactured solution

dofs = Matrix(zeros(nnodes*(dim+2),1))
testfunc = Matrix(zeros(nnodes*(dim+2),1))
for i in range(0,nnodes):
        for j in range(0,dim+2):
            dofs[i*(dim+2)+j] = U[i,j]
            testfunc[i*(dim+2)+j] = w[i,j]

## Compute LHS and RHS
print("\nCompute RHS\n")
rhs = Compute_RHS(rv.copy(), testfunc, do_simplifications)
rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)
    
print("\nCompute LHS\n")
lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS
lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)

## Residual for shock capturing
res_m_out = OutputMatrix_CollectingFactors(res_m, "res_m", mode)
res_e_out = OutputMatrix_CollectingFactors(res_e, "res_e", mode)

## Reading Template File
print("\nReading compressible_navier_stokes_cpp_template.cpp\n")
if(dim==2):
        templatefile = open("compressible_navier_stokes_cpp_template2D.cpp")
        outstring=templatefile.read()
        outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)
        outstring = outstring.replace("//substitute_res_m_2D", res_m_out)
        outstring = outstring.replace("//substitute_res_e_2D", res_e_out)
        ## Write the modified template
        print("\nWriting compressible_navier_stokes2D.cpp\n")
        out = open("compressible_navier_stokes2D.cpp",'w')
        out.write(outstring)
        out.close()
elif(dim == 3):
        templatefile = open("compressible_navier_stokes_cpp_template3D.cpp")
        outstring=templatefile.read()
        outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)
        outstring = outstring.replace("//substitute_res_m_3D", res_m_out)
        outstring = outstring.replace("//substitute_res_e_3D", res_e_out)

        ## Write the modified template
        print("\nWriting compressible_navier_stokes3D.cpp\n")
        out = open("compressible_navier_stokes3D.cpp",'w')
        out.write(outstring)
        out.close()
print("\nCompressible Navier Stokes Element Generated\n")