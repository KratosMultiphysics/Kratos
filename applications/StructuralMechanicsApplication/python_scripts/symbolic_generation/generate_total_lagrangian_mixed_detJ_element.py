import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Symbolic generation settings
mode = "c"
do_simplifications = False
dim = 2
nnodes = 3

if dim == 2:
    strain_size = 3
elif dim == 3:
    strain_size = 6
else:
    raise ValueError("Wrong dimension {}.".format(dim))
impose_partion_of_unity = False
N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

# Symbols definition
u = DefineMatrix('u',nnodes,dim) # Displacement (u(i,k) refers to the displacement of node i component k)
th = DefineVector('th',nnodes) # tetha variable representing the nodal det(J)
w = DefineMatrix('w',nnodes,dim) # Displacement test function
q = DefineVector('q',nnodes) # tetha test function
tau = sympy.Symbol("tau",positive=True) # Stabilization constant
rho0 = sympy.Symbol("rho0",positive=True) # Density in the initial configuration

S = DefineVector('S',strain_size) # Stress in Voigt notation (this will be returned by the constitutive law)
S = sympy.Matrix([[S[0],S[2]],[S[2],S[1]]]) # Definition of the stress tensor from the previous definition symbols

C = DefineSymmetricMatrix("C",strain_size,strain_size) # Constitutive matrix in Voigt notation (this will be returned by the constitutive law)
C = ConvertVoigtMatrixToTensor(C) # Definition of the 4th order constitutive tensor from the previous definition symbols

# Define the tetha interpolation at the Gauss point
th_gauss = 0
for n in range(nnodes):
  th_gauss += N[n]*th[n]

# Shape functions evaluation at the Gauss point
grad_w_gauss = w.transpose()*DN
grad_q_gauss = q.transpose()*DN
q_gauss = q.transpose()*N

# Define the deformation gradient tensor at the Gauss point
# F[i,j] = delta_ij + Du_i/Dx_j = delta_ij + sum_n DN_n*u_ni/Dx_j
F_gauss = sympy.Matrix(sympy.eye(dim,dim))
for i in range(dim):
  for j in range(dim):
    for n in range(nnodes):
        F_gauss[i,j] += DN[n,j]*u[n,i] #TODO: check that we are not using the transpose

# Define the Jacobian determinant at the Gauss point
j_gauss = sympy.det(F_gauss)

# Calculate the cofactor of the deformation gradient
# cof(F) = det(J)*F^{-T}
invF_gauss = F_gauss.inv()
cofF_gauss = j_gauss*(invF_gauss.transpose())






Fbar_gauss = 1/j_gauss**(1/3)*F_gauss
Fmod_gauss = th_gauss*Fbar_gauss

Emod_gauss = 1/2*(Fmod_gauss.transpose()*Fmod_gauss - sympy.eye(dim,dim))

grad_th_gauss = DN.transpose()*th

mom_first = DoubleContraction(grad_w_gauss, F_gauss* S)
mom_second = 0.0

mass_first=q_gauss*(j_gauss-th_gauss)

left = j_gauss*grad_q_gauss

tmp = DoubleContraction(C, (Fbar_gauss.transpose()*Fbar_gauss) ).tomatrix()
right = ( tau / (3.0*(th_gauss**(1/3))) ) * tmp * grad_th_gauss

mass_stab_term = left*right

functional = mom_first - mom_second+ mass_first[0] + mass_stab_term[0]