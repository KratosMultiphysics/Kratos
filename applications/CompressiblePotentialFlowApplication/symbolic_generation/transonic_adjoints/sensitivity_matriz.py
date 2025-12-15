##################################################################
#  Compressible potential-flow T3 element with density-upwind flux
#  and isentropic density law (symbolic)
#
#  Outputs symbolic expression for dR/dx (residual sensitivity wrt node coords)
#
#  Notes:
#   - Set use_isentropic=True to use the explicit isentropic rho(q2).
#   - rho0, a0 and gamma are symbolic reference quantities (you can substitute numeric values later).
##################################################################

import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Settings
mode = "c"
nnodes = 3
dim = 2
impose_partion_of_unity = False

# Choose density model
use_isentropic = True   # <-- set False to use a generic symbolic rho(q2)

# If using isentropic:
rho0 = sympy.symbols('rho0')   # reference density (symbolic)
a0   = sympy.symbols('a0')     # reference speed of sound (symbolic)
gamma = sympy.symbols('gamma') # ratio of specific heats

# Shape functions and local derivatives
N, DNDe = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

if dim == 2 and nnodes == 3:
    DNDe[0,0] = -1; DNDe[0,1] = -1
    DNDe[1,0] =  1; DNDe[1,1] =  0
    DNDe[2,0] =  0; DNDe[2,1] =  1

# Coordinates
x = DefineMatrix("x", nnodes, dim)

# Nodal potentials (interior)
p = DefineVector("p", nnodes)

# -- Geometry: Jacobian, DNDx, detJ
J     = sympy.simplify(x.transpose() * DNDe)
JInv  = sympy.simplify(J**(-1))
DNDx  = sympy.simplify(DNDe * JInv)
det_J = sympy.simplify(0.5 * J.det())

# -- Gradient of phi (interior)
grad_phi = sympy.Matrix([0,0])
for i in range(nnodes):
    grad_phi += p[i] * DNDx[i,:].transpose()
grad_phi = sympy.simplify(grad_phi)

# -- q^2 and rho definition ---------------------------------------
q2 = sympy.simplify(grad_phi[0]**2 + grad_phi[1]**2)

if use_isentropic:
    # Common isentropic density law (reference-based)
    # rho = rho0 * (1 + (gamma-1)/2 * q2 / a0^2)^{-1/(gamma-1)}
    rho_expr = sympy.simplify(rho0 * (1 + (gamma-1)/2 * q2 / (a0**2))**(-1/(gamma-1)))
    # we keep rho as an expression (so derivatives will expand properly)
    rho = sympy.simplify(rho_expr)
else:
    # Fully symbolic rho(q2) if you prefer not to fix the law here
    rho = sympy.Function("rho")(q2)

# -- Interior flux vector and its normal component (for an edge normal n)
F_int = sympy.simplify(rho * grad_phi)   # vector

# ---- Edges definition (node pairs) and helper functions -----------------------
edges = [(0,1), (1,2), (2,0)]

def edge_geometry(i, j):
    ex = x[j,0] - x[i,0]
    ey = x[j,1] - x[i,1]
    length = sympy.sqrt(ex**2 + ey**2)
    # choose a normal (element-local). For consistency pick [ey, -ex]
    nx = ey
    ny = -ex
    return sympy.simplify(ex), sympy.simplify(ey), sympy.simplify(length), sympy.simplify(nx), sympy.simplify(ny)

# ---- External states: for each edge define symbolic exterior grad phi and phi avg
gradPhiExt = [ sympy.Matrix([sympy.Function(f"gExt{e}_0")(), sympy.Function(f"gExt{e}_1")()]) for e in range(3) ]
phiExt    = [ sympy.Function(f"phiExt{e}")() for e in range(3) ]

# Interior average phi on each edge (linear -> average of nodal p on edge)
phiInt_edge = []
for (i,j) in edges:
    phiInt_edge.append( sympy.simplify( (p[i] + p[j]) / 2.0 ) )

# ---- Numerical flux (Lax-Friedrichs / Rusanov) on each edge -------------------
Fn_edges = []
for e_idx, (i,j) in enumerate(edges):
    ex, ey, edge_len, nx, ny = edge_geometry(i,j)
    nvec = sympy.Matrix([nx, ny])
    # interior normal flux component: Fi_n = rho_L * (grad_phi_L · n)
    Fi_n = sympy.simplify(F_int.dot(nvec))
    # exterior flux component (symbolic)
    gradR = gradPhiExt[e_idx]
    q2_R = sympy.simplify(gradR[0]**2 + gradR[1]**2)
    if use_isentropic:
        rhoR = sympy.simplify(rho0 * (1 + (gamma-1)/2 * q2_R / (a0**2))**(-1/(gamma-1)))
    else:
        rhoR = sympy.Function("rho")(q2_R)
    Fr_n = sympy.simplify(rhoR * (gradR.dot(nvec)))
    # alpha (wave speed approx) - use max(abs(Fi_n), abs(Fr_n)) (symbolic)
    alpha = sympy.simplify( sympy.Max(sympy.Abs(Fi_n), sympy.Abs(Fr_n)) )
    # numerical flux (scalar, normal component)
    Fn = sympy.simplify( 0.5*(Fi_n + Fr_n) - 0.5*alpha*( phiExt[e_idx] - phiInt_edge[e_idx] ) )
    Fn_edges.append( (Fn, edge_len, i, j) )

# ---- Volumetric residual (standard compressible divergence term) -------------
R_vol = sympy.Matrix(nnodes,1, lambda jj,_: sympy.simplify( DNDx[jj,:].dot(F_int) * det_J ))

# ---- Edge contributions to residual ------------------------------------------
R_edge = sympy.Matrix(nnodes, 1, lambda _,__ : 0)
for e_idx, (Fn, edge_len, ni, nj) in enumerate(Fn_edges):
    contrib = sympy.simplify(edge_len/2.0 * Fn)
    R_edge[ni,0] += contrib
    R_edge[nj,0] += contrib

# ---- Total residual vector R_total
R_total = sympy.simplify(R_vol + R_edge)

# ---- Sensitivity: derivative of R_total wrt node coordinates x(i,dim)
DB_Dx = DefineMatrix('DB_Dx', dim*nnodes, nnodes)
for inode in range(nnodes):
    for jdim in range(dim):
        for k in range(nnodes):
            DB_Dx[inode*dim + jdim, k] = sympy.simplify( sympy.diff( R_total[k,0], x[inode,jdim] ) )

# ---- Output (collect factors for C-like code)
dR_Dx_out = OutputMatrix_CollectingFactors(DB_Dx, "rOutput", mode)
print(dR_Dx_out)
