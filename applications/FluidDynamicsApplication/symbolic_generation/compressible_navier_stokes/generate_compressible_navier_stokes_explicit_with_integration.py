from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

from sympy import *
import pprint

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

# dim = params["dim"]         # Define Dimension in params.py
do_simplifications = False
mode = "c"                  # Output mode to a c++ file
is_explicit = True          # Explicit or implicit time integration

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

dim_vector = [2]
# dim_vector = [2,3]
subscales_vector = ["ASGS","OSS"]
# shock_capturing = "Isotropic"
shock_capturing = "Anisotropic"

for dim in dim_vector:
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
    r = DefineVector('r',n_nodes)                   # Sink term

    # Test functions defintiion
    w = DefineMatrix('w',n_nodes,block_size)	 # Variables field test

    # External terms definition
    f_ext = DefineMatrix('f_ext',n_nodes,dim)    # Forcing term #COMMENT for manufactured solution

    # Definition of other symbols
    if not is_explicit:
        # Backward differantiation coefficients
        bdf0 = Symbol('bdf0')
        bdf1 = Symbol('bdf1')
        bdf2 = Symbol('bdf2')
    k_sc = Symbol('k_sc', positive = True) # Shock capturing conductivity
    nu_sc = Symbol('nu_sc', positive = True) # Shock capturing viscosity
    k_st = Symbol('k_st', positive = True) # Stabilization conductivity approximation
    nu_st = Symbol('nu_st', positive = True) # Stabilization viscosity approximation
    lin_m = DefineVector('lin_m', dim) # Linearized momentum used in the anisotropic shock capturing matrices calculation
    lin_m_norm = Symbol('lin_m_norm', positive = True) # Linearized momentum norm. This has to be defined as a separated symbol to check it is non-zero.

    ### Construction of the variational equation
    Ug = DefineVector('Ug',block_size) # Dofs vector
    H = DefineMatrix('H',block_size,dim) # Gradient of U
    f = DefineVector('f',dim) # Body force vector
    rg = Symbol('rg', positive = True) # Source/Sink term
    V = DefineVector('V',block_size) # Test function
    Q = DefineMatrix('Q',block_size,dim) # Gradient of V
    acc = DefineVector('acc',block_size) # Derivative of Dofs/Time
    G = DefineMatrix('G',block_size,dim) # Diffusive Flux matrix
    res_proj = DefineVector('res_proj',block_size) # Residuals projection for the OSS

    ## Calculate the Gauss point residual
    ## Matrix Computation
    S = generate_source_term.computeS(f, rg, params)
    A = generate_convective_flux.computeA(Ug, params)
    if shock_capturing == "Isotropic":
        G = generate_diffusive_flux.ComputeDiffusiveFluxIsotropicShockCapturing(Ug, H, params, nu_sc, k_sc)
    elif shock_capturing == "Anisotropic":
        G = generate_diffusive_flux.ComputeDiffusiveFluxAnisotropicShockCapturing(Ug, H, params, nu_sc, k_sc, nu_st, k_st, lin_m, lin_m_norm)
    else:
        raise Exception("Wrong shock capturing method: \'" + shock_capturing + "\'. Available options are: \'Isotropic\' and \'Anisotropic\'")
    Tau = generate_stabilization_matrix.computeTau(params) # TODO: SOURCE TERMS CONSTANTS ARE NOT ADDED YET!

    ## Non-linear operator definition
    print("\nCompute Non-linear operator\n")
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
    # In the implicit case it is computed with the BDF formulas
    # In the explicit case it is linearised by using the values already stored in the database
    # In the explicit case it is required to add it to calculate the subscale intertial terms
    res = - acc - L

    # ## Isotropic Residual Based Shock Capturing
    # # Momentum residual
    # res_m = Matrix(zeros(dim,1))
    # for i in range(0,dim):
    #     res_m[i,0] = res[i+1,0]
    # # Energy residual
    # res_e = Matrix(zeros(1,1))
    # res_e[0,0] = res[dim+1]

    ## Non-linear adjoint operator definition
    print("\nCompute Non-linear Adjoint operator\n")
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
    diff_flux = 0.0
    n3 = Matrix(zeros(1,1))
    for j in range(dim):
        for k in range(block_size):
            diff_flux -= Q[k,j] * G[k,j]
    n3[0,0] = - diff_flux

    # Source term - FE scale
    n4 = V.transpose() * (S * Ug)

    # VMS_adjoint - Subscales
    subscales = DefineVector('subscales',block_size)
    # subscales = Tau * res
    n5 = L_adj.transpose() * subscales

    # Variational formulation (Galerkin functional)
    print("\nCompute Variational Formulation\n")
    if not is_explicit:
        rv = n1 + n2 + n3 + n4 + n5 # Implicit case (includes the inertial term n1)
    else:
        rv = n2 + n3 + n4 + n5 		# Explicit case

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
        if not is_explicit:
            # In the implicit case, calculate the time derivatives with the BDF2 formula
            acc_gauss = (bdf0 * U + bdf1 * Un + bdf2 * Unn).transpose()*N
        else:
            # In the explicit case, the acceleration is linearised taking the previous step one
            # Note that in the explicit case this acceleration is only used in the calculation of the stabilization terms
            acc_gauss = dUdt.transpose()*N
        r_gauss = (r.transpose()*N)[0]

        ## Gradients computation
        grad_U = DfjDxi(DN,U).transpose()

        ## Substitute the symbols in the residual
        SubstituteMatrixValue(res_gauss, Ug, U_gauss)
        SubstituteMatrixValue(res_gauss, acc, acc_gauss)
        SubstituteMatrixValue(res_gauss, H, grad_U)
        SubstituteMatrixValue(res_gauss, f, f_gauss)
        SubstituteScalarValue(res_gauss, rg, r_gauss)

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
        print("\nSubstitution of the discretized values at the gauss points")
        for i_gauss in range(n_gauss):
            print("\tGauss point: " + str(i_gauss))
            rv_gauss = rv.copy()

            ## Substitute the subscales model
            if subscales_type == "ASGS":
                asgs_subscales = Tau * res
                SubstituteMatrixValue(rv_gauss, subscales, asgs_subscales)
            elif subscales_type == "OSS":
                # oss_subscales = Tau * (res)
                oss_subscales = Tau * (res - res_proj)
                # oss_subscales = Tau * (res + res_proj)
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
            if not is_explicit:
                # In the implicit case, calculate the time derivatives with the BDF2 formula
                acc_gauss = (bdf0 * U + bdf1 * Un + bdf2 * Unn).transpose()*N
            else:
                # In the explicit case, the acceleration is linearised taking the previous step one
                # Note that in the explicit case this acceleration is only used in the calculation of the stabilization terms
                acc_gauss = dUdt.transpose()*N
            r_gauss = (r.transpose()*N)[0]

            ## Gauss pt. stabilization matrix calculation
            tau_gauss = generate_stabilization_matrix.computeTauOnGaussPoint(params, U_gauss)

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
        print("\nCompute RHS\n")
        rhs = Compute_RHS(rv_tot.copy(), testfunc, do_simplifications)
        rhs_out = OutputVector_CollectingFactors(rhs, "rRightHandSideBoundedVector", mode)

        if not is_explicit:
            print("\nCompute LHS\n")
            lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS
            lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)

        ## Reading and filling the template file
        print("\nSubstituting outstring in " + template_filename + " \n")
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
print("\nCompressible Navier Stokes Explicit Element Generated\n")