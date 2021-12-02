# import KratosMultiphysics
from sympy import Symbol, Matrix, zeros, diff

from KratosMultiphysics.sympy_fe_utilities import                  \
    DefineMatrix, DefineVector, DfjDxi, Compute_RHS, Compute_LHS,  \
    SubstituteMatrixValue, SubstituteScalarValue,                  \
    OutputMatrix_CollectingFactors, OutputVector_CollectingFactors

from params_dict import params
from shape_functions import DefineShapeFunctionsMatrix
import generate_convective_flux
import generate_diffusive_flux
import generate_source_term
import generate_stabilization_matrix

do_simplifications = False          # Simplify resulting differenctiations
subscales_vector = ["ASGS", "OSS"]  # Subscales types to be computed
dim_vector = [2]                    # Spatial dimensions to be computed
is_explicit = True                  # Explicit or implicit time integration
shock_capturing = True
mode = 'c'

template_filename = "" # TODO
with open(template_filename) as f:
    outstring = f.read()


def ComputeNonLinearOperator(block_size, Ug, S, A, H):
    L = Matrix(zeros(block_size, 1))
    for j in range(dim):
        # Convective operator product (A x grad(U))
        A_j = A[j]
        H_j = H.col(j)
        L += A_j * H_j
        # Diffusive flux
        # Note that the diffusive flux is not added as it will involve 2nd
        # order derivatives that vanish when introducing the linear FE
        # discretization

    # Source term addition
    L -= S * Ug

    return L


def ComputeAdjointOperator(dim, block_size, V, Q, S, A, H):
    L_adj = Matrix(zeros(block_size, 1))
    for j in range(dim):
        Q_j = Q.col(j)
        H_j = H.col(j)
        # Convective operator product
        A_j_trans = A[j].transpose()
        L_adj += A_j_trans * Q_j
        aux_conv = Matrix(zeros(block_size, block_size))
        for m in range(block_size):
            for n in range(block_size):
                A_j_trans_mn = A_j_trans[m, n]
                for k in range(block_size):
                    aux_conv[m, n] += diff(A_j_trans_mn, Ug[k]) * H_j[k]
        L_adj += aux_conv * V
        # Diffusive operator product
        # Note that the adjoint diffusive flux is not added as it will involve
        # 2nd order derivatives that vanish when introducing the linear FE
        # discretization

    # Source term addition
    L_adj += S.transpose() * V

    return L_adj


def ComputeVariationalFormulation(is_explicit, block_size, Q, G, A, H):
    # Convective term - FE scale
    conv_flux = zeros(block_size, 1)
    for j in range(dim):
        conv_flux += A[j] * H.col(j)
    n2 = - V.transpose() * conv_flux

    # Diffusive term - FE scale
    n3 = Matrix(zeros(1, 1))
    for j in range(dim):
        for k in range(block_size):
            n3[0, 0] += Q[k, j] * G[k, j]

    # Source term - FE scale
    n4 = V.transpose() * (S * Ug)

    # VMS_adjoint - Subscales
    subscales = DefineVector('subscales', block_size)
    n5 = L_adj.transpose() * subscales

    if is_explicit:
        rv = n2 + n3 + n4 + n5       # Explicit (without inertial term n1)
    else:
        n1 = - V.transpose()*acc  # Mass (inertial) term - FE scale
        rv = n1 + n2 + n3 + n4 + n5  # Implicit (includes the inertial term n1)
    return rv


def CalculateResidualsProjections(dim, n_nodes, n_gauss, mat_N, is_explicit):
    res_rho_proj = Matrix(zeros(n_nodes, 1))
    res_mom_proj = Matrix(zeros(n_nodes*dim, 1))
    res_tot_ener_proj = Matrix(zeros(n_nodes, 1))
    for i_gauss in range(n_gauss):
        print("\tGauss point: " + str(i_gauss))
        res_gauss = res.copy()

        # Get Gauss point geometry data
        N = mat_N.row(i_gauss)

        # Data interpolation at the gauss point
        U_gauss = N * U
        f_gauss = N * f_ext
        r_gauss = (N * r_ext)[0]
        mass_gauss = (N * m_ext)[0]
        if is_explicit:
            # previous step one. Note that in the explicit case this
            # acceleration is only used in the calculation of the stabilization
            # terms
            acc_gauss = N * dUdt
        else:
            # In the explicit case, the acceleration is linearised taking the
            # In the implicit case, calculate the time derivatives with the
            # BDF2 formula
            acc_gauss = N * (bdf0 * U + bdf1 * Un + bdf2 * Unn)

        # Gradients computation
        grad_U = DfjDxi(DN, U).transpose()

        # Substitute the symbols in the residual
        SubstituteMatrixValue(res_gauss, Ug, U_gauss)
        SubstituteMatrixValue(res_gauss, acc, acc_gauss)
        SubstituteMatrixValue(res_gauss, H, grad_U)
        SubstituteMatrixValue(res_gauss, f, f_gauss)
        SubstituteScalarValue(res_gauss, rg, r_gauss)
        SubstituteScalarValue(res_gauss, mg, mass_gauss)

        # Add the projection contributions
        for i_node in range(n_nodes):
            # Note that the weights will be added later on in the cpp
            res_rho_proj[i_node] += N[i_node] * res_gauss[0]
            res_tot_ener_proj[i_node] += N[i_node] * res_gauss[dim + 1]
            for d in range(dim):
                res_mom_proj[i_node * dim + d] += N[i_node] * res_gauss[1 + d]
    return (res_rho_proj, res_mom_proj, res_tot_ener_proj)


def OutputProjections(res_rho_proj, res_mom_proj, res_tot_ener_proj, outstring):
    res_rho_proj_out = OutputVector_CollectingFactors(res_rho_proj, "rho_proj", mode)
    res_mom_proj_out = OutputVector_CollectingFactors(res_mom_proj, "mom_proj", mode)
    res_tot_ener_proj_out = OutputVector_CollectingFactors(res_tot_ener_proj, "tot_ener_proj", mode)

    outstring = outstring.replace("//substitute_rho_proj_{}D".format(dim),
                                  res_rho_proj_out)
    outstring = outstring.replace("//substitute_mom_proj_{}D".format(dim),
                                  res_mom_proj_out)
    outstring = outstring.replace("//substitute_tot_ener_proj_{}D".format(dim),
                                  res_tot_ener_proj_out)
    return outstring


for dim in dim_vector:
    # Change dimension accordingly
    params["dim"] = dim

    # Shape functions and Gauss pts. settings
    (n_nodes, n_gauss) = {
        1: (1, 2),
        2: (4, 4),
        3: (8, 8)
    }[dim]

    DN = DefineMatrix('DN', n_nodes, dim)
    mat_N = DefineShapeFunctionsMatrix(dim, n_nodes, n_gauss)

    # Unknown fields definition (Used later for the gauss point interpolation)
    block_size = dim + 2
    U = DefineMatrix('U', n_nodes, block_size)  # Vector of Unknowns
    ResProj = DefineMatrix('ResProj', n_nodes, block_size)
    if is_explicit:  # Previous step data
        dUdt = DefineMatrix('dUdt', n_nodes, block_size)
    else:
        Un = DefineMatrix('Un', n_nodes, block_size)
        Unn = DefineMatrix('Unn', n_nodes, block_size)
        # Backward differantiation coefficients
        bdf0 = Symbol('bdf0')
        bdf1 = Symbol('bdf1')
        bdf2 = Symbol('bdf2')

    # Test functions defintiion
    w = DefineMatrix('w', n_nodes, block_size)   # Variables field test

    # External terms definition
    m_ext = DefineVector('m_ext', n_nodes)       # Mass source term
    r_ext = DefineVector('r_ext', n_nodes)       # Thermal sink/source term
    f_ext = DefineMatrix('f_ext', n_nodes, dim)  # Forcing term

    # Nodal artificial magnitudes
    alpha_sc_nodes = DefineVector('alpha_sc_nodes', n_nodes)  # mass diffusivity
    mu_sc_nodes = DefineVector('mu_sc_nodes', n_nodes)  # dynamic viscosity
    beta_sc_nodes = DefineVector('beta_sc_nodes', n_nodes)  # bulk viscosity
    lamb_sc_nodes = DefineVector('lamb_sc_nodes', n_nodes)  # bulk viscosity

    # Construction of the variational equation
    Ug = DefineVector('Ug', block_size)     # Dofs vector
    H = DefineMatrix('H', block_size, dim)  # Gradient of U
    mg = Symbol('mg')                       # Mass source term
    f = DefineVector('f', dim)              # Body force vector
    rg = Symbol('rg')                       # Thermal source/sink term
    V = DefineVector('V', block_size)       # Test function
    Q = DefineMatrix('Q', block_size, dim)  # Gradient of V
    acc = DefineVector('acc', block_size)   # Derivative of Dofs/Time
    G = DefineMatrix('G', block_size, dim)  # Diffusive Flux matrix
    res_proj = DefineVector('res_proj', block_size)  # OSS Residuals projection

    # Calculate the Gauss point residual
    # Matrix Computation
    S = generate_source_term.ComputeSourceMatrix(Ug, mg, f, rg, params)
    A = generate_convective_flux.ComputeEulerJacobianMatrix(Ug, params)
    # Shock capturing artificial magnitudes
    if shock_capturing:
        alpha_sc = Symbol('alpha_sc', positive=True)  # density diffusivity
        mu_sc = Symbol('mu_sc', positive=True)      # dynamic viscosity
        beta_sc = Symbol('beta_sc', positive=True)  # bulk viscosity
        lamb_sc = Symbol('lamb_sc', positive=True)  # thermal conductivity
        G = generate_diffusive_flux.ComputeDiffusiveFluxWithShockCapturing(
            Ug, H, params, alpha_sc, beta_sc, lamb_sc, mu_sc)
    else:
        G = generate_diffusive_flux.ComputeDiffusiveFlux(Ug, H, params)
    Tau = generate_stabilization_matrix.ComputeStabilizationMatrix(params)

    # Non-linear operator definition
    print("\nCompute non-linear operator\n")
    L = ComputeNonLinearOperator(block_size, Ug, S, A, H)

    # FE residuals definition
    # Note that we include the DOF time derivatives in both the implicit and
    # the explicit cases
    # It is required to include it in both cases to calculate the subscale
    # inertial component
    # In the implicit case it is computed with the BDF formulas
    # In the explicit case it is linearised by using the values already stored
    # in the database
    res = - acc - L

    # Non-linear adjoint operator definition
    print("\nCompute non-linear adjoint operator\n")
    L_adj = ComputeAdjointOperator(dim, block_size, V, Q, S, A, H)

    # Variational formulation (Galerkin functional)
    print("\nCompute variational formulation\n")
    rv = ComputeVariationalFormulation(is_explicit, block_size, Q, G, A, H)

    # OSS Residual projections calculation #
    # Calculate the residuals projection
    print("\nCalculate the projections of the residuals")
    projections = CalculateResidualsProjections(dim, n_nodes, n_gauss, mat_N, is_explicit)

    # Output the projections
    outstring = OutputProjections(*projections, outstring)

    # Algebraic form calculation #

