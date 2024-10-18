##################################################################
#    |  /           |
#    ' /   __| _` | __|  _ \   __|
#    . \  |   (   | |   (   |\__ `
#   _|\_\_|  \__,_|\__|\___/ ____/
#                   Multi-Physics
#
#  License:        BSD License
#                  Kratos default license: kratos/license.txt
#
#
#  Main authors:    Inigo Lopez and Riccardo Rossi

# This scripts computes the sensitivity matrix for the incompressible
# potential flow formulation simbolicaly using Kratos and sympy.
import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Symbolic generation settings
do_simplifications = False
# Output mode to a c++ file
mode = "c"
nnodes = 3
dim = 2

impose_partion_of_unity = False
# These are LOCAL derivatives
N, DNDe = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

#The node ordering corresponds with:
#     v
#     ^
#     |
#     2
#     |`\
#     |  `\
#     |    `\
#     |      `\
#     |        `\
#     0----------1 --> u

if(dim == 2 and nnodes==3):
    DNDe[0,0] = -1; DNDe[0,1] = -1
    DNDe[1,0] = 1; DNDe[1,1] = 0
    DNDe[2,0] = 0; DNDe[2,1] = 1

# Matrix with the coordinates  x(i,j) contains coordinates of node i component k
x = DefineMatrix('x', nnodes, dim)

# p[i] represents the phi of node i
p = DefineVector('p', nnodes)

# The Jacobian is calculated: derivative of the coordinates in terms of local derivatives
#  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
#J=|				|=	|			  |
#  |dy/dxi  dy/deta|	|y1-y0   y2-y0|
J = sympy.simplify(x.transpose()*DNDe)

# Inverse of the jacobian
JInverse = sympy.simplify((J)**(-1))

# Derivative of the shape functions w.r.t x and y
DNDx = sympy.simplify(DNDe*JInverse)

# Determinant of the jacobian (half of the area of the triangle)
det_J = sympy.simplify(0.5*J.det())

# Local stiffness matrix (before integrating over the area)
A = sympy.simplify(DNDx*DNDx.transpose())

# Multiplying with the area:
K = sympy.simplify(A*det_J)

# Computing the right hand side
B = sympy.simplify(-K * p)

# Deriving the right hand side w.r.t. each coordinate
DB_Dx = DefineMatrix('DB_Dx',dim*nnodes,nnodes)

# Loop over nodes: 1, 2, 3 for a triangle
for i_node in range(nnodes):
    # Loop over dimensions: x, y in 2d
    for j_dim in range(dim):
        # Looping over columns (i.e. )
        for k_x in range(nnodes):
            DB_Dx[ i_node * dim + j_dim , k_x ] = sympy.simplify(sympy.diff( B[k_x], x[i_node,j_dim]))

dB_Dx_out = OutputMatrix_CollectingFactors(DB_Dx, "rOutput", mode)

print(dB_Dx_out)