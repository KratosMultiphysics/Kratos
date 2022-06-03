import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Symbolic generation settings
mode = "c"
do_simplifications = False
dim = 2
n_nodes = dim + 1 # So far only simplex elements are considered
output_filename = "total_lagrangian_mixed_detJ_element.cpp"
template_filename = "total_lagrangian_mixed_detJ_element_cpp_template.cpp"

if dim == 2:
    strain_size = 3
elif dim == 3:
    strain_size = 6
else:
    raise ValueError("Wrong dimension {}.".format(dim))
impose_partion_of_unity = False
N,DN = DefineShapeFunctions(n_nodes, dim, impose_partion_of_unity)

# Symbols definition
u = DefineMatrix('u',n_nodes,dim) # Displacement (u(i,k) refers to the displacement of node i component k)
b = DefineMatrix('b',n_nodes,dim) # Body force (u(i,k) refers to the body force of node i component k)
th = DefineVector('th',n_nodes) # tetha variable representing the nodal det(J)
w = DefineMatrix('w',n_nodes,dim) # Displacement test function
q = DefineVector('q',n_nodes) # tetha test function
tau = sympy.Symbol("tau",positive=True) # Stabilization constant
rho0 = sympy.Symbol("rho0",positive=True) # Density in the initial configuration

S = DefineVector('S',strain_size) # Stress in Voigt notation (this will be returned by the constitutive law)
if dim == 2:
  S = sympy.Matrix([
    [S[0],S[2]],
    [S[2],S[1]]]) # Definition of the stress tensor from the previous definition symbols
else:
  S = sympy.Matrix([
    [S[0],S[3],S[5]],
    [S[3],S[1],S[4]],
    [S[5],S[4],S[2]]]) # Definition of the stress tensor from the previous definition symbols

C = DefineSymmetricMatrix("C",strain_size,strain_size) # Constitutive matrix in Voigt notation (this will be returned by the constitutive law)
C = ConvertVoigtMatrixToTensor(C) # Definition of the 4th order constitutive tensor from the previous definition symbols

# Define the tetha interpolations at the Gauss point
th_gauss = 0
for n in range(n_nodes):
  th_gauss += N[n]*th[n]
grad_th_gauss = DN.transpose()*th

# Define the body force interpolation at the Gauss point
b_gauss = b.transpose()*N

# Shape functions evaluation at the Gauss point
grad_w_gauss = w.transpose()*DN
grad_q_gauss = q.transpose()*DN
w_gauss = w.transpose()*N
q_gauss = q.transpose()*N

# Define the deformation gradient tensor at the Gauss point
# F[i,j] = delta_ij + Du_i/Dx_j = delta_ij + sum_n DN_n*u_ni/Dx_j
F_gauss = sympy.Matrix(sympy.eye(dim,dim))
for i in range(dim):
  for j in range(dim):
    for n in range(n_nodes):
        F_gauss[i,j] += DN[n,j]*u[n,i] #TODO: check that we are not using the transpose

# Define the Jacobian determinant at the Gauss point
j_gauss = sympy.det(F_gauss)

# Calculate the cofactor of the deformation gradient
# cof(F) = det(J)*F^{-T}
invF_gauss = F_gauss.inv()
cofF_gauss = j_gauss*(invF_gauss.transpose())

# Calculate the strain tensors
Fbar_gauss = (1/j_gauss**(1/dim))*F_gauss # Deviatoric deformation gradient tensor
Cbar_gauss = Fbar_gauss.transpose() * Fbar_gauss # Deviatoric right Cauchy-Green strain tensor
Ebar_gauss = 0.5*(Cbar_gauss - sympy.eye(dim,dim)) # Deviatoric Green strain tensor
# Fmod_gauss = th_gauss*Fbar_gauss
# Emod_gauss = 1/2*(Fmod_gauss.transpose()*Fmod_gauss - sympy.eye(dim,dim))

# Variational form
mom_first = DoubleContraction(grad_w_gauss, F_gauss* S)
mom_second = (w_gauss.transpose() * rho0 * b_gauss)[0]

mass_first = q_gauss[0] * (j_gauss - th_gauss)
tmp = (DoubleContraction(C, Cbar_gauss)).tomatrix() #FIXME: SHOULDN'T WE ONLY USE C IN THE DERIVATION OF THE LHS? HERE SHOULD BE THE STRESS.
aux_scalar = (tau / dim) * ((j_gauss / th_gauss)**(1.0/dim))
mass_stab_1 = (aux_scalar * grad_q_gauss * tmp * grad_th_gauss)[0]
mass_stab_2 = (tau * rho0 * grad_q_gauss * cofF_gauss.transpose() * b_gauss)[0]

functional = mom_first - mom_second + mass_first - (mass_stab_1 + mass_stab_2)
