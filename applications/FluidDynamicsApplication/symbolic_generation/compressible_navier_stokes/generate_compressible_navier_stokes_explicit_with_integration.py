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
    template_filename = "compressible_navier_stokes_cpp_template.cpp"
else:
    template_filename = "compressible_navier_stokes_explicit_cpp_template_with_integration.cpp"
templatefile = open(template_filename)
outstring = templatefile.read()

## Set the output filename
if template_filename == "compressible_navier_stokes_cpp_template":
    output_filename = "compressible_navier_stokes.cpp"
elif template_filename == "compressible_navier_stokes_explicit_cpp_template_with_integration.cpp":
    output_filename = "compressible_navier_stokes_explicit.cpp"
else:
    err_msg = "Wrong template_filename provided. Must be (template --> output):\n"
    err_msg +=  "\t- compressible_navier_stokes_cpp_template --> compressible_navier_stokes.cpp"
    err_msg +=  "\t- compressible_navier_stokes_explicit_cpp_template_with_integration.cpp --> compressible_navier_stokes_explicit.cpp"
    raise Exception(err_msg)

# dim_vector = [2]
dim_vector = [2,3]
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
    # Un = DefineMatrix('Un',n_nodes,block_size)      # Vector of Unknowns one step back (IMPLICIT ONLY)
    # Unn = DefineMatrix('Unn',n_nodes,block_size)    # Vector of Unknowns two steps back (IMPLICIT ONLY)
    r = DefineVector('r',n_nodes)                   # Sink term    #COMMENT for manufactured solution

    # Test functions defintiion
    w = DefineMatrix('w',n_nodes,block_size)	 # Variables field test

    # External terms definition
    f_ext = DefineMatrix('f_ext',n_nodes,dim)    # Forcing term #COMMENT for manufactured solution

    # Definition of other symbols
    # bdf0 = Symbol('bdf0')                       # Backward differantiation coefficients (IMPLICIT ONLY)
    # bdf1 = Symbol('bdf1')
    # bdf2 = Symbol('bdf2')
    v_sc = Symbol('v_sc')                       # Shock capturing Viscosity
    k_sc = Symbol('k_sc')                       # Shock capturing Conductivity

    ### Construction of the variational equation
    Ug = DefineVector('Ug',block_size)			# Dofs vector
    H = DefineMatrix('H',block_size,dim)		# Gradient of U
    f = DefineVector('f',dim)			        # Body force vector
    rg = Symbol('rg', positive = True)		    # Source/Sink term
    V = DefineVector('V',block_size)			# Test function
    Q = DefineMatrix('Q',block_size,dim)		# Gradient of V
    # acc = DefineVector('acc',block_size)        # Derivative of Dofs/Time (IMPLICIT ONLY)
    #G = DefineMatrix('G',block_size,dim)		# Diffusive Flux matrix
    Gsc = DefineMatrix('G',block_size,dim)      # Diffusive Flux matrix with Shock Capturing

    ## Calculate the Gauss point residual
    ## Matrix Computation
    S = generate_source_term.computeS(f, rg, params)
    #source_term.printS(S,params)
    A = generate_convective_flux.computeA(Ug, params)
    #convective_flux.printA(A,params)
    #G = diffusive_flux.computeG(Ug,params,H,G)
    Gsc = generate_diffusive_flux.computeGsc(Ug, params, H, Gsc, v_sc, k_sc)
    #diffusive_flux.printK(Gsc,params)
    Tau = generate_stabilization_matrix.computeTau(params)
    #stabilization_matrix.printTau(Tau, params)

    ## Nonlinear operator definition
    l1 = Matrix(zeros(dim+2,1))		            # Convective Matrix*Gradient of U
    A_small = []
    for j in range(0,dim):
        A_small = A[j]
        for ll in range(block_size):
            for mm in range(block_size):
                l1[ll] += A_small[ll,mm]*H[mm,j]

    l3 = S*Ug				                    # Source term
    print("\nCompute Non-linear operator\n")
    L = l1-l3                                   # Nonlinear operator

    ## Residual definition
    if not is_explicit:
        res = -acc - L # Implicit residual with inertial terms
    else:
        res = -L # Explicit residual without the inertial terms

    ## Nonlinear adjoint operator definition
    m1 = Matrix(zeros(dim+2,1))		            # Convective term
    psi = Matrix(zeros(dim+2,dim))

    for j in range(0,dim):
        A_T = A[j].transpose()
        for l in range(0,dim+2):
            for m in range(0,dim+2):
                psi[l,j] += A_T[l,m]*Q[m,j]
                for n in range(0,dim+2):
                    psi[l,j] +=diff(A_T[l,m],Ug[n])*H[n,j]*V[m]

    for s in range(0,dim+2):
        for j in range(0,dim):
            m1[s] += psi[s,j]

    m3 = S.transpose()*V			            # Source term

    L_adj = -m1-m3                              # Nonlinear adjoint operator

    ## Istotropic Residual Based Shock Capturing
    res_m = Matrix(zeros(dim,1))                # Momentum residual
    for i in range(0,dim):
        res_m[i,0] = res[i+1,0]

    res_e = Matrix(zeros(1,1))                  # Energy residual
    res_e[0,0] = res[dim+1]


    ## Variational Formulation - Final equation
    # n1 = V.transpose()*acc		                # Mass term - FE scale (IMPLICIT ONLY)

    temp = zeros(dim+2,1)
    A_smalll = []
    for i in range(0,dim):
        A_smalll = A[i]
        for ll in range(block_size):
            for mm in range(block_size):
                temp[ll] += A_smalll[ll,mm]*H[mm,i]

    n2 = V.transpose()*temp			            # Convective term - FE scale

    n3 = Matrix(zeros(1,1))                     # Diffusive term - FE scale

    for j in range(0,dim):
        for k in range(block_size):
            n3[0,0] += Q[k,j]*(-Gsc[k,j])       # G with shock capturing - FE scale

    n4 = -V.transpose()*(S*Ug)		            # Source term - FE scale

    n5 = L_adj.transpose()*(Tau*res)	        # VMS_adjoint - Subscales

    print("\nCompute Variational Formulation\n")
    # rv = n1+n2+n3+n4+n5 			            # VARIATIONAL FORMULATION - FINAL EQUATION (IMPLICIT)
    rv = n2 + n3 + n4 + n5 			            # VARIATIONAL FORMULATION - FINAL EQUATION (EXPLICIT)

    ### Substitution of the discretized values at the gauss points
    ## Loop and accumulate the residual in each Gauss point
    rv_tot = Matrix(zeros(1,1))
    res_m_tot = Matrix(zeros(dim,1))
    res_e_tot = Matrix(zeros(1,1))
    print("\nSubstitution of the discretized values at the gauss points")
    for i_gauss in range(n_gauss):
        print("\tGauss point: " + str(i_gauss))
        rv_gauss = rv.copy()
        res_m_gauss = res_m.copy()
        res_e_gauss = res_e.copy()

        ## Get Gauss point geometry data
        N = DefineVector('N', n_nodes)
        for i in range(n_nodes):
            N[i] = mat_N[i_gauss, i]

        ## Data interpolation at the gauss point
        U_gauss = U.transpose()*N
        w_gauss = w.transpose()*N
        f_gauss = f_ext.transpose()*N                     #COMMENT for manufactured solution
        # acc_gauss = (bdf0*U+bdf1*Un+bdf2*Unn).transpose()*N
        r_gauss = (r.transpose()*N)[0]                    #COMMENT for manufactured solution
        #r_gauss = Symbol('r_gauss', positive = True)     #USED fro manufactured solution

        ## Gauss pt. stabilization matrix calculation
        tau_gauss = generate_stabilization_matrix.computeTauOnGaussPoint(params, U_gauss)

        ## Gradients computation
        grad_U = DfjDxi(DN,U).transpose()
        grad_w = DfjDxi(DN,w).transpose()

        print("\t- Substitution in the variational formulation")
        SubstituteMatrixValue(rv_gauss, Ug, U_gauss)
        # SubstituteMatrixValue(rv_gauss, acc, acc_gauss) (IMPLICIT ONLY)
        SubstituteMatrixValue(rv_gauss, H, grad_U)
        SubstituteMatrixValue(rv_gauss, V, w_gauss)
        SubstituteMatrixValue(rv_gauss, Q, grad_w)
        SubstituteMatrixValue(rv_gauss, Tau, tau_gauss)
        SubstituteMatrixValue(rv_gauss, f, f_gauss)           #COMMENT for manufactured solution
        SubstituteScalarValue(rv_gauss, rg, r_gauss)          #COMMENT for manufactured solution

        # Note that the explicit shock capturing needs to be done in an external util
        # This is because we do not have (and do not want) to have the inertial terms involved in here
        if not is_explicit:
            print("\t- Substitution in the residual of momentum")
            SubstituteMatrixValue(res_m_gauss, Ug, U_gauss)
            # SubstituteMatrixValue(res_m_gauss, acc, acc_gauss) (IMPLICIT ONLY)
            SubstituteMatrixValue(res_m_gauss, H, grad_U)
            SubstituteMatrixValue(res_m_gauss, Tau, tau_gauss)
            SubstituteMatrixValue(res_m_gauss, f, f_gauss)       #COMMENT for manufactured solution
            SubstituteScalarValue(res_m_gauss, rg, r_gauss)      #COMMENT for manufactured solution

            print("\t- Substitution in the residual of total energy")
            SubstituteMatrixValue(res_e_gauss, Ug, U_gauss)
            # SubstituteMatrixValue(res_e_gauss, acc, acc_gauss) (IMPLICIT ONLY)
            SubstituteMatrixValue(res_e_gauss, H, grad_U)
            SubstituteMatrixValue(res_e_gauss, Tau, tau_gauss)
            SubstituteMatrixValue(res_e_gauss, f, f_gauss)       #COMMENT for manufactured solution
            SubstituteScalarValue(res_e_gauss, rg, r_gauss)      #COMMENT for manufactured solution

        ## Accumulate in the total value
        rv_tot += rv_gauss
        if not is_explicit:
            res_m_tot += res_m_gauss
            res_e_tot += res_e_gauss

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

    ## Residual for shock capturing
    if not is_explicit:
        res_m_out = OutputMatrix_CollectingFactors(res_m_tot, "res_m", mode)
        res_e_out = OutputMatrix_CollectingFactors(res_e_tot, "res_e", mode)

    ## Reading Template File
    print("\nReading compressible_navier_stokes_cpp_template.cpp\n")
    if(dim==2):
            outstring = outstring.replace("//substitute_rhs_2D", rhs_out)
            if not is_explicit:
                outstring = outstring.replace("//substitute_res_m_2D", res_m_out)
                outstring = outstring.replace("//substitute_res_e_2D", res_e_out)
                outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
    elif(dim == 3):
            outstring = outstring.replace("//substitute_rhs_3D", rhs_out)
            if not is_explicit:
                outstring = outstring.replace("//substitute_res_m_3D", res_m_out)
                outstring = outstring.replace("//substitute_res_e_3D", res_e_out)
                outstring = outstring.replace("//substitute_lhs_3D", lhs_out)

    ## Substitute the solution values container accesses
    for i_node in range(n_nodes):
        for j_block in range(block_size):
            to_substitute = 'U(' + str(i_node) + ',' + str(j_block) + ')'
            substituted_value = 'U_' + str(i_node) + '_' + str(j_block)
            outstring = outstring.replace(to_substitute, substituted_value)

    # ## Substitute the shape function container accesses
    # for i_gauss in range(n_gauss):
    #     for j_node in range(n_nodes):
    #         to_substitute = 'mat_N(' + str(i_gauss) + ',' + str(j_node) + ')'
    #         substituted_value = 'N_' + str(i_gauss) + '_' + str(j_node)
    #         outstring = outstring.replace(to_substitute, substituted_value)

    ## Substitute the shape function gradients container accesses
    for i_node in range(n_nodes):
        for j_dim in range(dim):
            to_substitute = 'DN(' + str(i_node) + ',' + str(j_dim) + ')'
            substituted_value = 'DN_DX_' + str(i_node) + '_' + str(j_dim)
            outstring = outstring.replace(to_substitute, substituted_value)

## Write the modified template
print("\nWriting compressible_navier_stokes_explicit.cpp\n")
out = open(output_filename,'w')
out.write(outstring)
out.close()
print("\nCompressible Navier Stokes Explicit Element Generated\n")