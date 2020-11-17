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
            mat_N[0,0] = 0.5
            mat_N[0,1] = 0.5
        elif n_gauss == 2:
            mat_N[0,0] = 0.788675
            mat_N[0,1] = 0.211325
            mat_N[1,0] = 0.211325
            mat_N[1,1] = 0.788675
        else:
            err_msg = "Invalid quadrature for dimension " + str(dim) + " and number of Gauss points " + str(n_gauss) + "."
    elif dim == 3:
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
    else:
        err_msg = "Invalid dimension " + str(dim) + "."
        raise Exception(err_msg)

    return mat_N

do_simplifications = False
mode = "c"                  # Output mode to a c++ file
is_explicit = True          # Explicit or implicit time integration

## Initialize the outstring to be filled with the template .cpp file
template_filename = "compressible_navier_stokes_wall_condition_cpp_template.cpp"
templatefile = open(template_filename)
outstring = templatefile.read()

## Set the output filename
if template_filename == "compressible_navier_stokes_wall_condition_cpp_template.cpp":
    output_filename = "compressible_navier_stokes_wall_condition.cpp"
else:
    err_msg = "Wrong template_filename provided. Must be (template --> output):\n"
    err_msg +=  "\t- compressible_navier_stokes_wall_condition_cpp_template.cpp --> compressible_navier_stokes_wall_condition.cpp"
    raise Exception(err_msg)

dim_vector = [2]
# dim_vector = [2,3]
# shock_capturing = "Isotropic"
# shock_capturing = "Anisotropic"
shock_capturing = "Peraire"

for dim in dim_vector:
    # Change dimension accordingly
    params["dim"] = dim

    # Shape functions and Gauss pts. settings
    if(dim == 2):
        n_nodes = 2
        n_gauss = 2
    elif(dim == 3):
        n_nodes = 3
        n_gauss = 3

    DN = DefineMatrix('DN', n_nodes, dim)
    mat_N = DefineShapeFunctionsMatrix(dim, n_nodes, n_gauss)

    # Unknown fields definition (Used later for the gauss point interpolation)
    block_size = dim + 2                            # Dimension of the vector of Unknowns
    U = DefineMatrix('U',n_nodes,block_size)	    # Vector of Unknowns (Density, Velocity[dim], Total Energy)

    # Test functions defintiion
    w = DefineMatrix('w',n_nodes,block_size)	 # Variables field test

    # Nodal artificial magnitudes
    beta_sc_nodes = DefineVector('beta_sc_nodes',n_nodes) # Nodal artificial bulk viscosity
    lamb_sc_nodes = DefineVector('lamb_sc_nodes',n_nodes) # Nodal artificial bulk viscosity

    # Definition of other symbols
    # k_sc = Symbol('k_sc', positive = True) # Shock capturing thermal diffusivity (TODO: Change to alpha)
    # nu_sc = Symbol('nu_sc', positive = True) # Shock capturing kinematic viscosity
    # k_st = Symbol('k_st', positive = True) # Stabilization thermal diffusivity approximation (TODO: Change to alpha)
    # nu_st = Symbol('nu_st', positive = True) # Stabilization kinematic viscosity approximation
    beta_sc = Symbol('beta_sc', positive = True) # Artificial bulk viscosity for shock capturing
    lamb_sc = Symbol('lamb_sc', positive = True) # Artificial thermal conductivity for shock capturing

    ### Construction of the variational equation
    Ug = DefineVector('Ug',block_size) # Dofs vector
    H = DefineMatrix('H',block_size,dim) # Gradient of U
    V = DefineVector('V',block_size) # Test function
    G = DefineMatrix('G',block_size,dim) # Diffusive Flux matrix

    ## Calculate the Gauss point diffusive flux
    # if shock_capturing == "Isotropic":
    #     G = generate_diffusive_flux.ComputeDiffusiveFluxIsotropicShockCapturing(Ug, H, params, nu_sc, k_sc)
    # elif shock_capturing == "Anisotropic":
    #     lin_m = DefineVector('lin_m', dim) # Linearized momentum used in the anisotropic shock capturing matrices calculation
    #     lin_m_norm = Symbol('lin_m_norm', positive = True) # Linearized momentum norm. This has to be defined as a separated symbol to check it is non-zero.
    #     G = generate_diffusive_flux.ComputeDiffusiveFluxAnisotropicShockCapturing(Ug, H, params, nu_sc, k_sc, nu_st, k_st, lin_m, lin_m_norm)
    # else:
    #     raise Exception("Wrong shock capturing method: \'" + shock_capturing + "\'. Available options are: \'Isotropic\' and \'Anisotropic\'")

    if shock_capturing == "Isotropic":
        G = generate_diffusive_flux.ComputeDiffusiveFluxIsotropicShockCapturing(Ug, H, params, nu_sc, k_sc)
    elif shock_capturing == "Peraire":
        G = generate_diffusive_flux.ComputeDiffusiveFluxPeraireShockCapturing(Ug, H, params, beta_sc, lamb_sc)
    elif shock_capturing == "Anisotropic":
        lin_m = DefineVector('lin_m', dim) # Linearized momentum used in the anisotropic shock capturing matrices calculation
        lin_m_norm = Symbol('lin_m_norm', positive = True) # Linearized momentum norm. This has to be defined as a separated symbol to check it is non-zero.
        G = generate_diffusive_flux.ComputeDiffusiveFluxAnisotropicShockCapturing(Ug, H, params, nu_sc, k_sc, nu_st, k_st, lin_m, lin_m_norm)
    else:
        raise Exception("Wrong shock capturing method: \'" + shock_capturing + "\'. Available options are: \'Peraire\', \'Isotropic\' and \'Anisotropic\'")

    # Variational formulation (boundary term)
    print("\nCompute variational boundary term\n")
    normal = DefineVector("normal",dim)
    ### TODO: CHECK IF THIS IS WETHER POSITVE OR NEGATIVE!!
    ### TODO: CHECK IF THIS IS WETHER POSITVE OR NEGATIVE!!
    ### TODO: CHECK IF THIS IS WETHER POSITVE OR NEGATIVE!!
    ### TODO: CHECK IF THIS IS WETHER POSITVE OR NEGATIVE!!
    ### TODO: CHECK IF THIS IS WETHER POSITVE OR NEGATIVE!!
    # rv = V.transpose() * (G * normal)
    rv = -V.transpose() * (G * normal)

    ### Substitution of the discretized values at the gauss points
    ## Loop and accumulate the residual in each Gauss point
    rv_tot = Matrix(zeros(1,1))

    print("\nSubstitution of the discretized values at the gauss points")
    for i_gauss in range(n_gauss):
        print("\tGauss point: " + str(i_gauss))
        rv_gauss = rv.copy()

        ## Get Gauss point geometry data
        N = DefineVector('N', n_nodes)
        for i in range(n_nodes):
            N[i] = mat_N[i_gauss, i]

        ## Data interpolation at the gauss point
        U_gauss = U.transpose() * N
        w_gauss = w.transpose() * N
        beta_sc_gauss = (beta_sc_nodes.transpose()*N)[0]
        lamb_sc_gauss = (lamb_sc_nodes.transpose()*N)[0]

        ## Gradients computation
        grad_U = DfjDxi(DN,U).transpose()
        grad_w = DfjDxi(DN,w).transpose()

        print("\t- Substitution in the variational formulation")
        SubstituteMatrixValue(rv_gauss, Ug, U_gauss)
        SubstituteMatrixValue(rv_gauss, H, grad_U)
        SubstituteMatrixValue(rv_gauss, V, w_gauss)
        SubstituteScalarValue(rv_gauss, beta_sc, beta_sc_gauss)
        SubstituteScalarValue(rv_gauss, lamb_sc, lamb_sc_gauss)

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
    outstring = outstring.replace("//substitute_rhs_" + str(dim) + "D", rhs_out)
    if not is_explicit:
        outstring = outstring.replace("//substitute_lhs_" + str(dim) + "D", lhs_out)

    # ## In the explicit element case the container values are referenced in the cpp to limit the container accesses to one per element
    # if is_explicit:
    #     ## Substitute the solution values container accesses
    #     for i_node in range(n_nodes):
    #         for j_block in range(block_size):
    #             to_substitute = 'U(' + str(i_node) + ',' + str(j_block) + ')'
    #             substituted_value = 'U_' + str(i_node) + '_' + str(j_block)
    #             outstring = outstring.replace(to_substitute, substituted_value)

    #     ## Substitute the solution values time derivatives container accesses
    #     for i_node in range(n_nodes):
    #         for j_block in range(block_size):
    #             to_substitute = 'dUdt(' + str(i_node) + ',' + str(j_block) + ')'
    #             substituted_value = 'dUdt_' + str(i_node) + '_' + str(j_block)
    #             outstring = outstring.replace(to_substitute, substituted_value)

    #     ## Substitute the shape function gradients container accesses
    #     for i_node in range(n_nodes):
    #         for j_dim in range(dim):
    #             to_substitute = 'DN(' + str(i_node) + ',' + str(j_dim) + ')'
    #             substituted_value = 'DN_DX_' + str(i_node) + '_' + str(j_dim)
    #             outstring = outstring.replace(to_substitute, substituted_value)

## Write the modified template
print("\nWriting " + output_filename + " \n")
out = open(output_filename,'w')
out.write(outstring)
out.close()
print("\nCompressible Navier-Stokes wall condition generated\n")