from sympy import *
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

from params_dict import params
import generate_convective_flux
import generate_diffusive_flux
import generate_source_term
import generate_stabilization_matrix

def DefineShapeFunctionsMatrix(dim, n_nodes, n_gauss):
    mat_N = DefineMatrix('mat_N', n_gauss, n_nodes)
    if dim == 2:
        if n_gauss == 1:
            mat_N[0,0] = 1.0 / 3.0
            mat_N[0,1] = 1.0 / 3.0
            mat_N[0,2] = 1.0 / 3.0
        elif n_gauss == 3:
            mat_N[0,0] = 2.0 / 3.0
            mat_N[0,1] = 1.0 / 6.0
            mat_N[0,2] = 1.0 / 6.0
            mat_N[1,0] = 1.0 / 6.0
            mat_N[1,1] = 2.0 / 3.0
            mat_N[1,2] = 1.0 / 6.0
            mat_N[2,0] = 1.0 / 6.0
            mat_N[2,1] = 1.0 / 6.0
            mat_N[2,2] = 2.0 / 3.0
        else:
            err_msg = "Invalid quadrature for dimension " + str(dim) + " and number of Gauss points " + str(n_gauss) + "."
    elif dim == 3:
        if n_gauss == 1:
            mat_N[0,0] = 1.0 / 4.0
            mat_N[0,1] = 1.0 / 4.0
            mat_N[0,2] = 1.0 / 4.0
            mat_N[0,3] = 1.0 / 4.0
        elif n_gauss == 4:
            mat_N[0,0] = 0.58541020
            mat_N[0,1] = 0.13819660
            mat_N[0,2] = 0.13819660
            mat_N[0,3] = 0.13819660
            mat_N[1,0] = 0.13819660
            mat_N[1,1] = 0.58541020
            mat_N[1,2] = 0.13819660
            mat_N[1,3] = 0.13819660
            mat_N[2,0] = 0.13819660
            mat_N[2,1] = 0.13819660
            mat_N[2,2] = 0.58541020
            mat_N[2,3] = 0.13819660
            mat_N[3,0] = 0.13819660
            mat_N[3,1] = 0.13819660
            mat_N[3,2] = 0.13819660
            mat_N[3,3] = 0.58541020
        else:
            err_msg = "Invalid quadrature for dimension " + str(dim) + " and number of Gauss points " + str(n_gauss) + "."
    else:
        err_msg = "Invalid dimension " + str(dim) + "."
        raise Exception(err_msg)

    return mat_N

mode = "c"                          # Output mode to a c++ file
is_explicit = True                  # Explicit (True) or implicit (False) time integration
do_simplifications = False          # Simplify resulting differenctiations
dim_vector = [2,3]                  # Spatial dimensions to be computed
shock_capturing = True              # Add physics-based shock capturing contribution
subscales_vector = ["ASGS","OSS"]   # Subscales types to be computed

## Initialize the outstring to be filled with the template .cpp file
if not is_explicit:
    template_filename = "compressible_navier_stokes_cpp_template_with_integration.cpp"
else:
    template_filename = "compressible_navier_stokes_explicit_cpp_template_with_integration.cpp"
templatefile = open(template_filename)
outstring = templatefile.read()

## Set the output filename
if template_filename == "compressible_navier_stokes_cpp_template_with_integration.cpp":
    output_filename = "compressible_navier_stokes.cpp"
elif template_filename == "compressible_navier_stokes_explicit_cpp_template_with_integration.cpp":
    output_filename = "compressible_navier_stokes_explicit.cpp"
else:
    err_msg = "Wrong template_filename provided. Must be (template --> output):\n"
    err_msg +=  "\t- compressible_navier_stokes_cpp_template_with_integration.cpp --> compressible_navier_stokes.cpp"
    err_msg +=  "\t- compressible_navier_stokes_explicit_cpp_template_with_integration.cpp --> compressible_navier_stokes_explicit.cpp"
    raise Exception(err_msg)

for dim in dim_vector:
    print("\nComputing dimension: " + str(dim) + "D\n")

    # Change dimension accordingly
    params["dim"] = dim

    # Shape functions and Gauss pts. settings
    if(dim == 2):
        n_nodes = 3
        n_gauss = 3
    elif(dim == 3):
        n_nodes = 4
        n_gauss = 4

    DN = DefineMatrix('DN', n_nodes, dim)
    mat_N = DefineShapeFunctionsMatrix(dim, n_nodes, n_gauss)

    # Unknown fields definition (Used later for the gauss point interpolation)
    block_size = dim + 2                            # Dimension of the vector of Unknowns
    U = DefineMatrix('U',n_nodes,block_size)	    # Vector of Unknowns (Density, Velocity[dim], Total Energy)
    ResProj = DefineMatrix('ResProj',n_nodes,block_size)	# Vector of residuals projection
    if not is_explicit:
        Un = DefineMatrix('Un',n_nodes,block_size)      # Vector of Unknowns one step back
        Unn = DefineMatrix('Unn',n_nodes,block_size)    # Vector of Unknowns two steps back
    else:
        dUdt = DefineMatrix('dUdt',n_nodes,block_size)  # Vector of Unknowns time derivatives (Density, Velocity[dim], Total Energy)

    # Test functions defintiion
    w = DefineMatrix('w',n_nodes,block_size)	 # Variables field test

    # External terms definition
    m_ext = DefineVector('m_ext',n_nodes)        # Mass source term
    r_ext = DefineVector('r_ext',n_nodes)        # Thermal sink/source term
    f_ext = DefineMatrix('f_ext',n_nodes,dim)    # Forcing term

    # Nodal artificial magnitudes
    mu_sc_nodes = DefineVector('mu_sc_nodes',n_nodes) # Nodal artificial dynamic viscosity
    beta_sc_nodes = DefineVector('beta_sc_nodes',n_nodes) # Nodal artificial bulk viscosity
    lamb_sc_nodes = DefineVector('lamb_sc_nodes',n_nodes) # Nodal artificial bulk viscosity

    # Definition of other symbols
    if not is_explicit:
        # Backward differantiation coefficients
        bdf0 = Symbol('bdf0')
        bdf1 = Symbol('bdf1')
        bdf2 = Symbol('bdf2')

    ### Construction of the variational equation
    Ug = DefineVector('Ug',block_size) # Dofs vector
    H = DefineMatrix('H',block_size,dim) # Gradient of U
    mg = Symbol('mg') # Mass source term
    f = DefineVector('f',dim) # Body force vector
    rg = Symbol('rg') # Thermal source/sink term
    V = DefineVector('V',block_size) # Test function
    Q = DefineMatrix('Q',block_size,dim) # Gradient of V
    acc = DefineVector('acc',block_size) # Derivative of Dofs/Time
    G = DefineMatrix('G',block_size,dim) # Diffusive Flux matrix
    res_proj = DefineVector('res_proj',block_size) # Residuals projection for the OSS

    ## Calculate the Gauss point residual
    ## Matrix Computation
    S = generate_source_term.ComputeSourceMatrix(Ug, mg, f, rg, params)
    A = generate_convective_flux.ComputeEulerJacobianMatrix(Ug, params)
    if shock_capturing:
        mu_sc = Symbol('mu_sc', positive = True) # Artificial dynamic viscosity for shock capturing
        beta_sc = Symbol('beta_sc', positive = True) # Artificial bulk viscosity for shock capturing
        lamb_sc = Symbol('lamb_sc', positive = True) # Artificial thermal conductivity for shock capturing
        G = generate_diffusive_flux.ComputeDiffusiveFluxWithPhysicsBasedShockCapturing(Ug, H, params, beta_sc, lamb_sc, mu_sc)
    else:
        G = generate_diffusive_flux.ComputeDiffusiveFlux(Ug, H, params)
    Tau = generate_stabilization_matrix.ComputeStabilizationMatrix(params)

    ## Non-linear operator definition
    print("\nCompute non-linear operator\n")
    L = Matrix(zeros(block_size,1))
    for j in range(dim):
        # Convective operator product (A x grad(U))
        A_j = A[j]
        H_j = H.col(j)
        L += A_j * H_j
        # Diffusive flux
        # Note that the diffusive flux is not added as it will involve 2nd order derivatives that vanish when introducing the linear FE discretization
    # Source term addition
    L -= S * Ug

    ## FE residuals definition
    # Note that we include the DOF time derivatives in both the implicit and the explicit cases
    # It is required to include it in both cases to calculate the subscale inertial component
    # In the implicit case it is computed with the BDF formulas
    # In the explicit case it is linearised by using the values already stored in the database
    res = - acc - L

    ## Non-linear adjoint operator definition
    print("\nCompute non-linear adjoint operator\n")
    L_adj = Matrix(zeros(block_size,1))
    for j in range(dim):
        Q_j = Q.col(j)
        H_j = H.col(j)
        # Convective operator product
        A_j_trans = A[j].transpose()
        L_adj += A_j_trans * Q_j
        aux_conv = Matrix(zeros(block_size, block_size))
        for m in range(block_size):
            for n in range(block_size):
                A_j_trans_mn = A_j_trans[m,n]
                for l in range(block_size):
                    aux_conv[m,n] += diff(A_j_trans_mn, Ug[l]) * H_j[l]
        L_adj += aux_conv * V
        # Diffusive operator product
        # Note that the adjoint diffusive flux is not added as it will involve 2nd order derivatives that vanish when introducing the linear FE discretization

    # Source term addition
    L_adj += S.transpose() * V

    ## Variational Formulation - Final equation
    # Mass (inertial) term - FE scale (only computed in the implicit case)
    if not is_explicit:
        n1 = - V.transpose()*acc

    # Convective term - FE scale
    conv_flux = zeros(block_size, 1)
    for j in range(dim):
        conv_flux += A[j] * H.col(j)
    n2 = - V.transpose() * conv_flux

    # Diffusive term - FE scale
    n3 = Matrix(zeros(1,1))
    for j in range(dim):
        for k in range(block_size):
            n3[0,0] += Q[k,j] * G[k,j]

    # Source term - FE scale
    n4 = V.transpose() * (S * Ug)

    # VMS_adjoint - Subscales
    subscales = DefineVector('subscales',block_size)
    n5 = L_adj.transpose() * subscales

    # Variational formulation (Galerkin functional)
    print("\nCompute variational formulation\n")
    if not is_explicit:
        rv = n1 + n2 + n3 + n4 + n5 # Implicit case (includes the inertial term n1)
    else:
        rv = n2 + n3 + n4 + n5 		# Explicit case (without inertial term n1)

    #### OSS Residual projections calculation ####
    # Calculate the residuals projection
    print("\nCalculate the projections of the residuals")
    res_rho_proj = Matrix(zeros(n_nodes,1))
    res_mom_proj = Matrix(zeros(n_nodes*dim,1))
    res_tot_ener_proj = Matrix(zeros(n_nodes,1))
    for i_gauss in range(n_gauss):
        print("\tGauss point: " + str(i_gauss))
        res_gauss = res.copy()

        ## Get Gauss point geometry data
        N = DefineVector('N', n_nodes)
        for i_node in range(n_nodes):
            N[i_node] = mat_N[i_gauss, i_node]

        ## Data interpolation at the gauss point
        U_gauss = U.transpose() * N
        f_gauss = f_ext.transpose() * N
        r_gauss = (r_ext.transpose()*N)[0]
        mass_gauss = (m_ext.transpose()*N)[0]
        if not is_explicit:
            # In the implicit case, calculate the time derivatives with the BDF2 formula
            acc_gauss = (bdf0 * U + bdf1 * Un + bdf2 * Unn).transpose()*N
        else:
            # In the explicit case, the acceleration is linearised taking the previous step one
            # Note that in the explicit case this acceleration is only used in the calculation of the stabilization terms
            acc_gauss = dUdt.transpose()*N

        ## Gradients computation
        grad_U = DfjDxi(DN,U).transpose()

        ## Substitute the symbols in the residual
        SubstituteMatrixValue(res_gauss, Ug, U_gauss)
        SubstituteMatrixValue(res_gauss, acc, acc_gauss)
        SubstituteMatrixValue(res_gauss, H, grad_U)
        SubstituteMatrixValue(res_gauss, f, f_gauss)
        SubstituteScalarValue(res_gauss, rg, r_gauss)
        SubstituteScalarValue(res_gauss, mg, mass_gauss)

        ## Add the projection contributions
        for i_node in range(n_nodes):
            # Note that the weights will be added later on in the cpp
            res_rho_proj[i_node] += N[i_node] * res_gauss[0]
            for d in range(dim):
                res_mom_proj[i_node * dim + d] += N[i_node] * res_gauss[1 + d]
            res_tot_ener_proj[i_node] += N[i_node] * res_gauss[dim + 1]

    ## Output the projections
    res_rho_proj_out = OutputVector_CollectingFactors(res_rho_proj, "rho_proj", mode)
    res_mom_proj_out = OutputVector_CollectingFactors(res_mom_proj, "mom_proj", mode)
    res_tot_ener_proj_out = OutputVector_CollectingFactors(res_tot_ener_proj, "tot_ener_proj", mode)
    outstring = outstring.replace("//substitute_rho_proj_"+ str(dim) +"D", res_rho_proj_out)
    outstring = outstring.replace("//substitute_mom_proj_"+ str(dim) +"D", res_mom_proj_out)
    outstring = outstring.replace("//substitute_tot_ener_proj_"+ str(dim) +"D", res_tot_ener_proj_out)

    #### Algebraic form calculation ####
    for subscales_type in subscales_vector:
        ### Substitution of the discretized values at the gauss points
        ## Loop and accumulate the residual in each Gauss point
        rv_tot = Matrix(zeros(1,1))

        print("\nSubscales type: " + subscales_type)
        print("\n- Substitution of the discretized values at the gauss points")
        for i_gauss in range(n_gauss):
            print("\tGauss point: " + str(i_gauss))
            rv_gauss = rv.copy()

            ## Substitute the subscales model
            if subscales_type == "ASGS":
                asgs_subscales = Tau * res
                SubstituteMatrixValue(rv_gauss, subscales, asgs_subscales)
            elif subscales_type == "OSS":
                oss_subscales = Tau * (res - res_proj)
                SubstituteMatrixValue(rv_gauss, subscales, oss_subscales)
            else:
                raise Exception("Wrong subscales type!")

            ## Get Gauss point geometry data
            N = DefineVector('N', n_nodes)
            for i in range(n_nodes):
                N[i] = mat_N[i_gauss, i]

            ## Data interpolation at the gauss point
            U_gauss = U.transpose() * N
            w_gauss = w.transpose() * N
            f_gauss = f_ext.transpose() * N
            r_gauss = (r_ext.transpose()*N)[0]
            mass_gauss = (m_ext.transpose()*N)[0]
            mu_sc_gauss = (mu_sc_nodes.transpose()*N)[0]
            beta_sc_gauss = (beta_sc_nodes.transpose()*N)[0]
            lamb_sc_gauss = (lamb_sc_nodes.transpose()*N)[0]
            if not is_explicit:
                # In the implicit case, calculate the time derivatives with the BDF2 formula
                acc_gauss = (bdf0 * U + bdf1 * Un + bdf2 * Unn).transpose()*N
            else:
                # In the explicit case, the acceleration is linearised taking the previous step one
                # Note that in the explicit case this acceleration is only used in the calculation of the stabilization terms
                acc_gauss = dUdt.transpose()*N

            ## Gauss pt. stabilization matrix calculation
            if shock_capturing:
                tau_gauss = generate_stabilization_matrix.ComputeStabilizationMatrixOnGaussPoint(params, U_gauss, f_gauss, r_gauss, mu_sc_gauss, lamb_sc_gauss)
            else:
                tau_gauss = generate_stabilization_matrix.ComputeStabilizationMatrixOnGaussPoint(params, U_gauss, f_gauss, r_gauss)

            ## If OSS, residual projections interpolation
            if subscales_type == "OSS":
                res_proj_gauss = ResProj.transpose() * N

            ## Gradients computation
            grad_U = DfjDxi(DN,U).transpose()
            grad_w = DfjDxi(DN,w).transpose()

            print("\t- Substitution in the variational formulation")
            SubstituteMatrixValue(rv_gauss, Ug, U_gauss)
            SubstituteMatrixValue(rv_gauss, acc, acc_gauss)
            SubstituteMatrixValue(rv_gauss, H, grad_U)
            SubstituteMatrixValue(rv_gauss, V, w_gauss)
            SubstituteMatrixValue(rv_gauss, Q, grad_w)
            SubstituteMatrixValue(rv_gauss, Tau, tau_gauss)
            SubstituteMatrixValue(rv_gauss, f, f_gauss)
            SubstituteScalarValue(rv_gauss, rg, r_gauss)
            SubstituteScalarValue(rv_gauss, mg, mass_gauss)
            SubstituteScalarValue(rv_gauss, mu_sc, mu_sc_gauss)
            SubstituteScalarValue(rv_gauss, beta_sc, beta_sc_gauss)
            SubstituteScalarValue(rv_gauss, lamb_sc, lamb_sc_gauss)
            if subscales_type == "OSS":
                SubstituteMatrixValue(rv_gauss, res_proj, res_proj_gauss)

            ## Accumulate in the total value
            rv_tot += rv_gauss

        ## Set the DOFs and test function matrices to do the differentiation
        dofs = Matrix(zeros(n_nodes*(dim+2),1))
        testfunc = Matrix(zeros(n_nodes*(dim+2),1))
        for i in range(0,n_nodes):
                for j in range(0,dim+2):
                    dofs[i*(dim+2)+j] = U[i,j]
                    testfunc[i*(dim+2)+j] = w[i,j]

        ## Compute LHS and RHS
        print("\n- Compute RHS")
        rhs = Compute_RHS(rv_tot.copy(), testfunc, do_simplifications)
        rhs_out = OutputVector_CollectingFactors(rhs, "rRightHandSideBoundedVector", mode)

        if not is_explicit:
            print("\n- Compute LHS")
            lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS
            lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)

        ## Reading and filling the template file
        print("\n- Substituting outstring in " + template_filename + " \n")
        outstring = outstring.replace("//substitute_rhs_" + str(dim) + "D_" + subscales_type, rhs_out)
        if not is_explicit:
            outstring = outstring.replace("//substitute_lhs_" + str(dim) + "D_" + subscales_type, lhs_out)

        ## In the explicit element case the container values are referenced in the cpp to limit the container accesses to one per element
        if is_explicit:
            ## Substitute the solution values container accesses
            for i_node in range(n_nodes):
                for j_block in range(block_size):
                    to_substitute = 'U(' + str(i_node) + ',' + str(j_block) + ')'
                    substituted_value = 'U_' + str(i_node) + '_' + str(j_block)
                    outstring = outstring.replace(to_substitute, substituted_value)

            ## Substitute the solution values time derivatives container accesses
            for i_node in range(n_nodes):
                for j_block in range(block_size):
                    to_substitute = 'dUdt(' + str(i_node) + ',' + str(j_block) + ')'
                    substituted_value = 'dUdt_' + str(i_node) + '_' + str(j_block)
                    outstring = outstring.replace(to_substitute, substituted_value)

            ## Substitute the shape function gradients container accesses
            for i_node in range(n_nodes):
                for j_dim in range(dim):
                    to_substitute = 'DN(' + str(i_node) + ',' + str(j_dim) + ')'
                    substituted_value = 'DN_DX_' + str(i_node) + '_' + str(j_dim)
                    outstring = outstring.replace(to_substitute, substituted_value)

            ## Substitute the residuals projection container accesses
            for i_node in range(n_nodes):
                for j_dim in range(block_size):
                    to_substitute = 'ResProj(' + str(i_node) + ',' + str(j_dim) + ')'
                    substituted_value = 'ResProj_' + str(i_node) + '_' + str(j_dim)
                    outstring = outstring.replace(to_substitute, substituted_value)

## Write the modified template
print("\nWriting " + output_filename + " \n")
out = open(output_filename,'w')
out.write(outstring)
out.close()
print("\n" + output_filename + " generated\n")