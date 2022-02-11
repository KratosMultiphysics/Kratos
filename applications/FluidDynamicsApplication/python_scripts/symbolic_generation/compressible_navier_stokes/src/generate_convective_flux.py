import sympy
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes \
    .src.quantity_converter import QuantityConverter


def ComputeEulerJacobianMatrix(dofs, params, primitives):
    """This function calculates the Euler Jacobian matrix for convection"""

    # Auxiliary variables
    dim = params.dim

    rho = QuantityConverter.density(primitives, params)
    e_tot = QuantityConverter.total_energy(primitives, params, rho=rho)

    vel = primitives.V
    p = primitives.P

    # Define and fill the convective flux matrix E
    E = sympy.zeros(dim + 2, dim, real=True)
    delta = sympy.eye(dim)
    E[0, :]    = rho * vel.T
    E[1:-1, :] = rho*vel*vel.T + p*delta
    E[-1, :]   = vel.T * (e_tot + p)

    # Intermediate step to obtain dE/dV
    dE_dV = []
    primitives_vector = primitives.PrimitivesVector()
    for j in range(dim):
        Astar_j = sympy.zeros(dim + 2, dim + 2)
        for m in range(dim + 2):
            for n in range(dim + 2):
                Astar_j[m, n] = sympy.diff(E[m, j], primitives_vector[n])
        dE_dV.append(Astar_j)

    # Obtain the Euler Jacobian Matrix A := dE/dU = dE/dV * dV/dU
    dV_dU = primitives.dVdU(dofs)
    A = dE_dV * dV_dU
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
