##################################################################
#  Symbolic compressible potential-flow sensitivity matrix (T3)
#
#  Author: ChatGPT (following Kratos style and your previous script)
#
##################################################################

import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Symbolic generation settings
do_simplifications = False
mode = "c"          # output mode
nnodes = 3
dim = 2

impose_partion_of_unity = False
N, DNDe = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

# --- TRIANGLE LOCAL COORDINATES -------------------------------------------------
if dim == 2 and nnodes == 3:
    DNDe[0,0] = -1; DNDe[0,1] = -1
    DNDe[1,0] =  1; DNDe[1,1] =  0
    DNDe[2,0] =  0; DNDe[2,1] =  1

# Coordinates (symbolic)
x = DefineMatrix("x", nnodes, dim)

# DOFs: potential phi at nodes
p = DefineVector("p", nnodes)

# ---- COMPUTE GEOMETRIC QUANTITIES --------------------------------------------
J      = sympy.simplify(x.transpose() * DNDe)
JInv   = sympy.simplify(J ** (-1))
DNDx   = sympy.simplify(DNDe * JInv)
det_J  = sympy.simplify(0.5 * J.det())

# ---- GRADIENT OF POTENTIAL ----------------------------------------------------
# grad(phi) = sum_i p[i] * grad(N_i)
grad_phi = sympy.Matrix([0,0])
for i in range(nnodes):
    grad_phi += p[i] * DNDx[i,:].transpose()

# ---- SYMBOLIC COMPRESSIBLE DENSITY FUNCTION ----------------------------------
# Here density = rho(q^2), with q^2 = |grad(phi)|^2
q2 = sympy.simplify(grad_phi[0]**2 + grad_phi[1]**2)

# Leave rho generic (symbolic function)
rho = sympy.Function("rho")(q2)

# ---- FLUX VECTOR --------------------------------------------------------------
#   F = rho * grad(phi)
F = sympy.simplify(rho * grad_phi)

# ---- RESIDUAL (strong form projected to FE basis) -----------------------------
# Element residual: R_j = ∫ ( ∇N_j · F ) dA
# For T3, integrand is constant → R = B * det(J)
R = sympy.Matrix(nnodes, 1, lambda j,_ : DNDx[j,:] * F)   # ∇N_j · F
R = sympy.simplify(R * det_J)

# ---- SENSITIVITY: derivative wrt coordinates ----------------------------------
dR_Dx = DefineMatrix("dR_Dx", dim*nnodes, nnodes)

for i_node in range(nnodes):
    for j_dim in range(dim):
        for k in range(nnodes):
            dR_Dx[i_node*dim + j_dim, k] = sympy.simplify(
                sympy.diff(R[k], x[i_node, j_dim])
            )

# ---- OUTPUT CODE --------------------------------------------------------------
out = OutputMatrix_CollectingFactors(dR_Dx, "rOutput", mode)
print(out)
