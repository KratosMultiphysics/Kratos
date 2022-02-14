import sympy
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .src.quantity_converter import QuantityConverter


def ComputeEulerJacobianMatrix(U, primitives, params):
    """This function calculates the Euler Jacobian matrix for convection"""

    # Auxiliary variables
    dim = params.dim
    blocksize = dim + 2

    rho = QuantityConverter.density(primitives, params)
    e_tot = QuantityConverter.total_energy(primitives, params, rho=rho)

    vel = primitives.V
    p = primitives.P

    # Define and fill the convective flux matrix E
    E = sympy.zeros(blocksize, dim, real=True)
    delta = sympy.eye(dim)
    E[0, :]    = rho * vel.T
    E[1:-1, :] = rho*vel*vel.T + p*delta
    E[-1, :]   = vel.T * (e_tot + p)

    # Obtain the Euler Jacobian Matrix A := dE/dU = dE/dV * dV/dU
    V = primitives.AsVector()
    dV_dU = QuantityConverter.dVdU(U, primitives, params)
    A = [0] * dim
    for j in range(dim):
        dE_dV = sympy.zeros(blocksize, blocksize)

        for i in range(blocksize):
            for l in range(blocksize):
                dE_dV[i, l] = sympy.diff(E[i, j], V[l])
        A[j] = dE_dV * dV_dU[j, :]

    A.simplify()
    return A


def printA(A, params):
    dim = params.dim
    tmp = []
    print("The convective matrix is:\n")
    for j in range(0, dim):
        tmp = A[j]
        for i in range(0, dim + 2):
            for k in range(0, dim + 2):
                print("A[{},{},{}]=[{}]\n".format(j, i, k, tmp[i, k]))

    return 0
