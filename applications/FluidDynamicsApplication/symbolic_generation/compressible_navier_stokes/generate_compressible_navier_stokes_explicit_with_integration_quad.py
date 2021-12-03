# import KratosMultiphysics
from sympy import var, Symbol, Matrix, zeros, diff

from KratosMultiphysics.sympy_fe_utilities import                  \
    DfjDxi, Compute_RHS, Compute_LHS,  \
    SubstituteMatrixValue, SubstituteScalarValue,                  \
    OutputMatrix_CollectingFactors, OutputVector_CollectingFactors

import generate_convective_flux
import generate_diffusive_flux
import generate_source_term
import generate_stabilization_matrix


#################################################
#                                               #
#                 DATA ENTRY                    #
#                                               #
#################################################


do_simplifications = False          # Simplify resulting differenctiations
subscales_vector = ["ASGS", "OSS"]  # Subscales types to be computed
dim_vector = [2]                    # Spatial dimensions to be computed
is_explicit = True                  # Explicit or implicit time integration
shock_capturing = True
mode = 'c'

template_filename = "compressible_navier_stokes_explicit_cpp_quad_template_with_integration.cpp"
output_filename = "compressible_navier_stokes_quad.cpp"


#################################################
#                                               #
#              SUPPORT FUNCTIONS                #
#                                               #
#################################################


def DefineMatrix(name, n, m):
    return Matrix(n, m, lambda i, j: Symbol("{}({},{})".format(name, i, j), real=True))


def DefineVector(name, n):
    return Matrix(n, 1, lambda i, _: Symbol("{}({})".format(name, i), real=True))


def ZeroMatrix(rows, cols):
    return Matrix(rows, cols, lambda *_: 0.0)


def ZeroVector(rows):
    return ZeroMatrix(rows, 1)


def ComputeNonLinearOperator():
    L = ZeroVector(block_size)
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


def ComputeAdjointOperator():
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


def ComputeVariationalFormulation():
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
    return (rv, subscales)


def CalculateResidualsProjections():
    res_rho_proj = ZeroVector(n_nodes)
    res_mom_proj = ZeroVector(n_nodes*dim)
    res_tot_ener_proj = ZeroVector(n_nodes)

    res_gauss = res.copy()

    # Data interpolation at the gauss point
    U_gauss = U.T * N
    f_gauss = f_ext.T * N
    r_gauss = (r_ext.T * N)[0]
    mass_gauss = (m_ext.T * N)[0]
    if is_explicit:
        # previous step one. Note that in the explicit case this
        # acceleration is only used in the calculation of the stabilization
        # terms
        acc_gauss = dUdt.T * N
    else:
        # In the explicit case, the acceleration is linearised taking the
        # In the implicit case, calculate the time derivatives with the
        # BDF formula
        acc_gauss = (bdf[0] * U + bdf[1] * Un + bdf[2] * Unn).T * N

    # Gradients computation
    grad_U = DfjDxi(DN_DX, U).transpose()

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
    res_rho_proj = OutputVector_CollectingFactors(res_rho_proj, "rho_proj", mode, replace_indices=False, assignment_op=" += ")
    res_mom_proj = OutputVector_CollectingFactors(res_mom_proj, "mom_proj", mode, replace_indices=False, assignment_op=" += ")
    res_ener_proj = OutputVector_CollectingFactors(res_tot_ener_proj, "tot_ener_proj", mode, replace_indices=False, assignment_op=" += ")

    outstring = outstring.replace("//substitute_rho_proj_{}D".format(dim),
                                  res_rho_proj)
    outstring = outstring.replace("//substitute_mom_proj_{}D".format(dim),
                                  res_mom_proj)
    outstring = outstring.replace("//substitute_tot_ener_proj_{}D".format(dim),
                                  res_ener_proj)
    return outstring


def ComputeLeftAndRightHandSide():
    # Substitution of the discretized values at the gauss points
    # Loop and accumulate the residual in each Gauss point
    print("\nSubscales type: " + subscales_type)
    print("\n- Substitution of the discretized values at the gauss points")

    rv_ = rv.copy()  # Copy for this subscales type

    # Substitute the subscales model
    if subscales_type == "ASGS":
        asgs_subscales = Tau * res
        SubstituteMatrixValue(rv_, subscales, asgs_subscales)
    elif subscales_type == "OSS":
        oss_subscales = Tau * (res - res_proj)
        SubstituteMatrixValue(rv_, subscales, oss_subscales)
    else:
        raise ValueError("Unknown subscales type: {}".format(subscales_type))

    # Data interpolation at the gauss point
    U_gauss = U.transpose() * N
    w_gauss = w.transpose() * N
    f_gauss = f_ext.transpose() * N
    r_gauss = (r_ext.transpose()*N)[0]
    mass_gauss = (m_ext.transpose()*N)[0]
    alpha_sc_gauss = (alpha_sc_nodes.transpose()*N)[0]
    mu_sc_gauss = (mu_sc_nodes.transpose()*N)[0]
    beta_sc_gauss = (beta_sc_nodes.transpose()*N)[0]
    lamb_sc_gauss = (lamb_sc_nodes.transpose()*N)[0]

    if is_explicit:
        # In the explicit case, the acceleration is linearised taking the
        # previous step one. Note that in the explicit case this acceleration
        # is only used in the calculation of the stabilization terms
        acc_gauss = dUdt.transpose()*N
    else:
        # In the implicit case, calculate the time derivatives with the BDF2
        # formula
        acc_gauss = (bdf[0] * U + bdf[1] * Un + bdf[2] * Unn).transpose()*N

    # Gauss pt. stabilization matrix calculation
    if shock_capturing:
        tau_gauss = generate_stabilization_matrix.ComputeStabilizationMatrixOnGaussPoint(params, U_gauss, f_gauss, r_gauss, mu_sc_gauss, lamb_sc_gauss)
    else:
        tau_gauss = generate_stabilization_matrix.ComputeStabilizationMatrixOnGaussPoint(params, U_gauss, f_gauss, r_gauss)

    # If OSS, residual projections interpolation
    if subscales_type == "OSS":
        res_proj_gauss = ResProj.transpose() * N

    # Gradients computation
    grad_U = DfjDxi(DN_DX, U).transpose()
    grad_w = DfjDxi(DN_DX, w).transpose()

    print("\t- Substitution in the variational formulation")
    SubstituteMatrixValue(rv_, Ug, U_gauss)
    SubstituteMatrixValue(rv_, acc, acc_gauss)
    SubstituteMatrixValue(rv_, H, grad_U)
    SubstituteMatrixValue(rv_, V, w_gauss)
    SubstituteMatrixValue(rv_, Q, grad_w)
    SubstituteMatrixValue(rv_, Tau, tau_gauss)
    SubstituteMatrixValue(rv_, f, f_gauss)
    SubstituteScalarValue(rv_, rg, r_gauss)
    SubstituteScalarValue(rv_, mg, mass_gauss)
    SubstituteScalarValue(rv_, alpha_sc, alpha_sc_gauss)
    SubstituteScalarValue(rv_, mu_sc, mu_sc_gauss)
    SubstituteScalarValue(rv_, beta_sc, beta_sc_gauss)
    SubstituteScalarValue(rv_, lamb_sc, lamb_sc_gauss)
    if subscales_type == "OSS":
        SubstituteMatrixValue(rv_, res_proj, res_proj_gauss)

    # Compute LHS and RHS

    # Set the DOFs and test function matrices to do the differentiation
    dofs = ZeroMatrix(n_nodes*block_size, 1)
    testfunc = ZeroMatrix(n_nodes*block_size, 1)
    for i in range(n_nodes):
        for j in range(block_size):
            dofs[i*block_size + j] = U[i, j]
            testfunc[i*block_size + j] = w[i, j]

    print("\n- Compute RHS")
    rhs = Compute_RHS(rv_, testfunc, do_simplifications)

    lhs = None
    if not is_explicit:
        print("\n- Compute LHS")
        lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications)

    return lhs, rhs


def WriteLeftAndRightHandSide(lhs, rhs, outstring):
    rhs_out = OutputVector_CollectingFactors(rhs, "rRightHandSideBoundedVector", mode, replace_indices=False, assignment_op=" += ")

    if not is_explicit:
        lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode, replace_indices=False, assignment_op=" += ")

    # Reading and filling the template file
    print("\n- Substituting outstring in " + template_filename + " \n")
    outstring = outstring.replace("//substitute_rhs_" + str(dim) + "D_" + subscales_type, rhs_out)
    if not is_explicit:
        outstring = outstring.replace("//substitute_lhs_" + str(dim) + "D_" + subscales_type, lhs_out)
    return outstring


#################################################
#                                               #
#                     MAIN                      #
#                                               #
#################################################


with open(template_filename) as f:
    outstring = f.read()

params = {
    "dim": -1,                                       # Dimension
    "mu": Symbol('data.mu', positive=True),          # Dynamic viscosity
    "h": Symbol('data.h', positive=True),            # Element size
    "lambda": Symbol('data.lambda', positive=True),  # Thermal Conductivity
    "c_v": Symbol('data.c_v', positive=True),        # Specific Heat at constant volume
    "gamma": Symbol('data.gamma', positive=True),    # Gamma (Cp/Cv)
    "stab_c1": Symbol('stab_c1', positive=True),     # Algorithm constant
    "stab_c2": Symbol('stab_c2', positive=True),     # Algorithm constant
    "stab_c3": Symbol('stab_c3', positive=True),     # Algorithm constant
}


for dim in dim_vector:
    # Change dimension accordingly
    params["dim"] = dim

    # Shape functions and Gauss pts. settings
    (n_nodes, n_gauss) = {
        1: (1, 2),  # Line
        2: (4, 4),  # Quad
        3: (8, 8)   # Hexa
    }[dim]

    N = DefineVector("N", n_nodes)
    DN_DX = DefineMatrix("DN_DX", n_nodes, dim)

    # Unknown fields definition (Used later for the gauss point interpolation)
    block_size = dim + 2
    U = DefineMatrix('data.U', n_nodes, block_size)  # Vector of Unknowns
    ResProj = DefineMatrix('data.ResProj', n_nodes, block_size)
    if is_explicit:  # Previous step data
        dUdt = DefineMatrix('data.dUdt', n_nodes, block_size)
    else:
        Un = DefineMatrix('data.Un', n_nodes, block_size)
        Unn = DefineMatrix('data.Unn', n_nodes, block_size)
        # Backward differantiation coefficients
        bdf = [Symbol('bdf0'), Symbol('bdf1'), Symbol('bdf2')]

    # Test functions defintiion
    w = DefineMatrix('w', n_nodes, block_size)   # Variables field test

    # External terms definition
    m_ext = DefineVector('data.m_ext', n_nodes)       # Mass source term
    r_ext = DefineVector('data.r_ext', n_nodes)       # Thermal sink/source term
    f_ext = DefineMatrix('data.f_ext', n_nodes, dim)  # Forcing term

    # Nodal artificial magnitudes
    alpha_sc_nodes = DefineVector('data.alpha_sc_nodes', n_nodes)  # mass diffusivity
    mu_sc_nodes = DefineVector('data.mu_sc_nodes', n_nodes)        # dynamic viscosity
    beta_sc_nodes = DefineVector('data.beta_sc_nodes', n_nodes)    # bulk viscosity
    lamb_sc_nodes = DefineVector('data.lamb_sc_nodes', n_nodes)    # bulk viscosity

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
    L = ComputeNonLinearOperator()

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
    L_adj = ComputeAdjointOperator()

    # Variational formulation (Galerkin functional)
    print("\nCompute variational formulation\n")
    (rv, subscales) = ComputeVariationalFormulation()

    # OSS Residual projections calculation #
    # Calculate the residuals projection
    print("\nCalculate the projections of the residuals")
    projections = CalculateResidualsProjections()
    outstring = OutputProjections(*projections, outstring)

    # LSH and RHS
    for subscales_type in subscales_vector:
        matrices = ComputeLeftAndRightHandSide()
        outstring = WriteLeftAndRightHandSide(*matrices, outstring)


print("\nWriting " + output_filename + " \n")
with open(output_filename, 'w') as f:
    f.write(outstring)
print("\n" + output_filename + " generated\n")