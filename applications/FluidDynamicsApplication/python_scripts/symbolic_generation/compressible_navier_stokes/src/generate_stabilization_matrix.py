import sympy
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src.defines \
    import CompressibleNavierStokesDefines as defs


def ComputeStabilizationMatrix(params):
    """This function calculates the stabilization matrix"""

    dim = params.dim                # Spatial dimensions
    Tau = defs.ZeroMatrix(dim + 2, dim + 2)   # Stabilization Matrix

    tau1 = sympy.Symbol('tau1')
    tau2 = sympy.Symbol('tau2')
    tau3 = sympy.Symbol('tau3')

    Tau[0, 0] = tau1
    for i in range(0,dim):
        Tau[i + 1, i + 1] = tau2
    Tau[dim + 1, dim + 1] = tau3
    return(Tau)


def ComputeStabilizationMatrixOnGaussPoint(params, U_gauss, f_gauss, r_gauss, mu_sc_gauss=0.0, lamb_sc_gauss=0.0):
    """This function calculates the stabilization matrix on a Gauss point"""

    # Calculate auxiliary values
    rho_g = U_gauss[0]
    e_t_g = U_gauss[params.dim + 1]
    norm_v_squared = 0.0
    norm_f_squared = 0.0
    for d in range(params.dim):
        norm_v_squared += (U_gauss[d + 1] * U_gauss[d + 1]) / (rho_g * rho_g)
        norm_f_squared += f_gauss[d] * f_gauss[d]
    norm_v = sympy.sqrt(norm_v_squared)
    nu = (params.mu + mu_sc_gauss) / rho_g
    alpha = (params.lamb + lamb_sc_gauss) / (rho_g * params.gamma * params.c_v)

    # Calculate sound speed
    c = sympy.sqrt(params.gamma * (params.gamma - 1) * ((e_t_g / rho_g) - ((1.0 / 2.0) * norm_v_squared)))

    # Calculate stabilization constants
    tau1_inv = (params.stab_c2 * (norm_v + c)) / params.h
    tau1_inv += params.stab_c3 * sympy.sqrt((r_gauss**2 + 2.0 * c**2 * norm_f_squared + sympy.sqrt(r_gauss**4 + 4.0 * c**2 * norm_f_squared * r_gauss**2)) / (2.0 * c**4))
    tau2_inv = ((params.stab_c1 * 4.0 * nu) / (3 * params.h**2)) + tau1_inv
    tau3_inv = (params.stab_c1 * alpha / params.h**2) + tau1_inv

    # Save the obtained values in the stabilization matrix
    Tau = defs.ZeroMatrix(params.dim + 2, params.dim + 2)
    Tau[0, 0] = 1.0 / tau1_inv
    for i in range(params.dim):
        Tau[i + 1, i + 1] = 1.0 / tau2_inv
    Tau[params.dim + 1, params.dim + 1] = 1.0 / tau3_inv

    return(Tau)


def PrintStabilizationMatrix(Tau, params):
    """Auxiliary function to print the stabilization matrix"""

    print("The Stabilization term matrix is:\n")
    dim = params.dim
    for i in range(0, dim + 2):
        for j in range(0, dim + 2):
            print("Tau[", i, ",", j, "]=", Tau[i, j], "\n")

    return 0
