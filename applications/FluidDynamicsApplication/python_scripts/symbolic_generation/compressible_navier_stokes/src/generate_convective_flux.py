import sympy
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .src.defines import CompressibleNavierStokesDefines as defs


def ComputeEulerJacobianMatrix(dofs, params, primitives):
    """This function calculates the Euler Jacobian matrix for convection"""

    # Auxiliary variables
    dim = params.dim
    mom = dofs[1:-1]
    e_tot = dofs[-1]
    vel = primitives.V
    p = primitives.P

    # Define and fill the convective flux matrix
    E = defs.Matrix('E', dim + 2, dim, real=True)
    for j in range(dim):
        E[0, j] = mom[j]
        for d in range(dim):
            E[1 + d, j] = mom[j]*vel[d]
            if j == d:
                E[1 + d, j] += p
        E[dim + 1, j] = vel[j] * (e_tot + p)

    # Define and derive the Euler Jacobian matrix
    A = []
    for j in range(dim):
        A_j = defs.Matrix('A_j', dim + 2, dim + 2, real=True)
        for m in range(dim + 2):
            for n in range(dim + 2):
                A_j[m, n] = sympy.diff(E[m, j], dofs[n])
        A.append(A_j)

    return A


def printA(A, params):
    dim = params.dim
    tmp = []
    print("The convective matrix is:\n")
    for j in range(0, dim):
        tmp = A[j]
        for i in range(0, dim + 2):
            for k in range(0, dim + 2):
                print("A[{},{},{}]=[]\n".format(j, i, k, tmp[i, k]))

    return 0
