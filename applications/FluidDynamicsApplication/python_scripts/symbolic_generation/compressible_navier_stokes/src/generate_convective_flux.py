import sympy
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .src.defines import CompressibleNavierStokesDefines as defs


def ComputeEulerJacobianMatrix(dofs, params):
    """This function calculates the Euler Jacobian matrix for convection"""

    # Auxiliary variables
    dim = params.dim
    gamma = params.gamma
    rho_0 = params.rho_0
    A_JWL = params.A_JWL
    B_JWL = params.B_JWL
    omega = params.omega
    R1 = params.R1
    R2 = params.R2
    rho = dofs[0]
    mom = []
    vel = []
    mom_prod = 0.0
    for i in range(dim):
        mom.append(dofs[i + 1])
        vel.append(dofs[i + 1] / rho)
        mom_prod += dofs[i + 1]**2
    e_tot = dofs[dim + 1]
    V = rho_0 / rho  # Relative volume
    e_internal = e_tot - 0.5 * mom_prod / rho # Interna Energy

    p = A_JWL * (1 - omega / (R1 * V)) * sympy.exp(-R1 * V) + \
    B_JWL * (1 - omega / (R2 * V)) * sympy.exp(-R2 * V) + \
    (omega * e_internal) / V


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
