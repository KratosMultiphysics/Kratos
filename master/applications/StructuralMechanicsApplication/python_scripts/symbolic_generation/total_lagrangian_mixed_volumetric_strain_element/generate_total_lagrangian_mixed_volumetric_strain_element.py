import sympy
from KratosMultiphysics import *
from KratosMultiphysics import kratos_utilities
from KratosMultiphysics.sympy_fe_utilities import *

# Symbolic generation settings
mode = "c"
dim_vect = [2, 3]
do_simplifications = False
output_filename = "total_lagrangian_mixed_volumetric_strain_element.cpp"
template_filename = "total_lagrangian_mixed_volumetric_strain_element_template.cpp"
tokens_filenames = []

for dim in dim_vect:
    n_nodes = dim + 1  # So far only simplex elements are considered
    block_size = dim + 1
    local_size = n_nodes * block_size

    if dim == 2:
        strain_size = 3
    elif dim == 3:
        strain_size = 6
    else:
        raise ValueError("Wrong dimension {}.".format(dim))

    impose_partion_of_unity = False
    N,DN = DefineShapeFunctions(n_nodes, dim, impose_partion_of_unity)

    # Symbols definition
    u = DefineMatrix('u',n_nodes,dim) # Displacement (u(i,k) refers to the displacement of node i component k)
    b = DefineMatrix('b',n_nodes,dim) # Body force (u(i,k) refers to the body force of node i component k)
    th = DefineVector('th',n_nodes) # Variable representing the nodal det(J) - 1 (volumetric deformation increment)
    w = DefineMatrix('w',n_nodes,dim) # Displacement test function
    q = DefineVector('q',n_nodes) # Deformation Jacobian test function
    w_g = sympy.Symbol("w_g", positive=True) # Gauss point integration weight
    tau_u = sympy.Symbol("tau_u",positive=True) # Displacement stabilization constant
    tau_th = sympy.Symbol("tau_th",positive=True) # Deformation Jacobian stabilization constant

    S = DefineVector('S',strain_size) # Stress in Voigt notation (this will be returned by the constitutive law)
    if dim == 2:
        S = sympy.Matrix([
            [S[0],S[2]],
            [S[2],S[1]]]) # Definition of the stress tensor from the previous definition symbols
    else:
        S = sympy.Matrix([
            [S[0],S[3],S[5]],
            [S[3],S[1],S[4]],
            [S[5],S[4],S[2]]]) # Definition of the stress tensor from the previous definition symbols

    C = DefineSymmetricMatrix("C",strain_size,strain_size) # Constitutive matrix in Voigt notation (this will be returned by the constitutive law)
    C = ConvertVoigtMatrixToTensor(C) # Definition of the 4th order constitutive tensor from the previous definition symbols

    # Define the tetha interpolations at the Gauss point
    th_gauss = 0
    for n in range(n_nodes):
        th_gauss += N[n]*th[n]
        grad_th_gauss = DN.transpose()*th

    # Define the body force interpolation at the Gauss point
    b_gauss = DefineVector('b_gauss',dim)

    # Shape functions evaluation at the Gauss point
    grad_w_gauss = w.transpose()*DN
    grad_q_gauss = q.transpose()*DN
    w_gauss = w.transpose()*N
    q_gauss = q.transpose()*N

    # Define the deformation gradient tensor at the Gauss point
    # Note that this is defined only in terms of the displacement (not enriched). See article Remark 1
    # F[i,j] = delta_ij + Du_i/Dx_j = delta_ij + sum_n DN_n*u_ni/Dx_j
    F_gauss = sympy.Matrix(sympy.eye(dim,dim))
    for i in range(dim):
        for j in range(dim):
            for n in range(n_nodes):
                F_gauss[i,j] += DN[n,j]*u[n,i]

    # Define the Jacobian determinant at the Gauss point
    j_gauss = sympy.det(F_gauss)

    # Calculate the cofactor of the deformation gradient
    # cof(F) = det(J)*F^{-T}
    invF_gauss = F_gauss.inv(method="LU")
    cofF_gauss = j_gauss*(invF_gauss.transpose())

    # Calculate the strain tensors
    # Note that for the multiplicative decomposition of the strain we assume the dimension to be the volumetric deformation exponent
    # As a consequence, current implementation is valid for the general 3D case and the 2D plain strain one (not valid for the 2D plane stress)
    Fbar_gauss = (1/j_gauss**sympy.Rational(1,dim))*F_gauss # Deviatoric deformation gradient tensor
    Cbar_gauss = Fbar_gauss.transpose() * Fbar_gauss # Deviatoric right Cauchy-Green strain tensor
    Ebar_gauss = 0.5*(Cbar_gauss - sympy.eye(dim,dim)) # Deviatoric Green-Lagrange strain tensor

    # Calculate the equivalent strain tensors
    Fmod_gauss = ((1.0 + th_gauss)**sympy.Rational(1,dim))*Fbar_gauss # Equivalent (enriched) deformation gradient tensor
    Cmod_gauss = Fmod_gauss.transpose() * Fmod_gauss # Equivalent (enriched) right Cauchy-Green strain tensor
    Emod_gauss = 0.5*(Cmod_gauss - sympy.eye(dim,dim)) # Equivalent (enriched) Green-Lagrange strain tensor

    # Variational form
    tmp = (DoubleContraction(C, F_gauss.transpose()*F_gauss)).tomatrix()

    mom_first = DoubleContraction(grad_w_gauss, F_gauss* S)
    mom_second = (w_gauss.transpose() * b_gauss)[0]
    mom_aux_scalar = (tau_th / dim) * ((1+th_gauss)**sympy.Rational(2-dim,dim)) * (1.0/j_gauss**sympy.Rational(2,dim))
    mom_stab = DoubleContraction(grad_w_gauss, mom_aux_scalar * (1 + th_gauss - j_gauss) * tmp)

    mass_first = (1.0 - tau_th) * q_gauss[0] * (1.0+th_gauss - j_gauss)
    mass_aux_scalar = (tau_u / dim) * ((j_gauss / (1.0+th_gauss))**sympy.Rational(dim-2,dim))
    mass_stab_1 = (mass_aux_scalar * grad_q_gauss * tmp * grad_th_gauss)[0]
    mass_stab_2 = (tau_u * grad_q_gauss * cofF_gauss.transpose() * b_gauss)[0]

    functional = mom_second - mom_first + mom_stab + mass_first + mass_stab_1 + mass_stab_2
    functional_array = sympy.Matrix([functional])

    # Define DOFs and test function vectors
    dofs = sympy.zeros(local_size, 1)
    testfunc = sympy.zeros(local_size, 1)
    for i in range(n_nodes):
        # Displacement DOFs and test functions
        for k in range(dim):
            dofs[i*block_size + k] = u[i,k]
            testfunc[i*block_size + k] = w[i,k]
        # Jacobian determinant DOFs and test functions
        dofs[i*block_size + dim] = th[i,0]
        testfunc[i*block_size + dim] = q[i,0]

    # Compute RHS (functional differentiation w.r.t. the shape functions)
    # Note that the stress is included as a symbolic variable, which is assumed to computed by the constitutive law module
    rhs = Compute_RHS(functional_array.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(w_g*rhs, "rRightHandSideVector", mode, indentation_level=2, assignment_op="+=")

    # Compute LHS (RHS(residual) differentiation w.r.t. the DOFs)
    # Note that 'S' (stress symbolic variable) is substituted by a definition of a function in terms of E for the LHS differentiation
    # Otherwise the displacement terms inside the modified Green strain would not be considered in the differentiation
    # Also note that a direct substitution by C:E is not valid in this case as this would imply a wrong differentiation of the LHS terms involving S (e.g. geometric stiffness)

    # Create an auxiliary symbol for the Green-Lagrange strain as a function such that E(u,tetha)
    E = sympy.MatrixSymbol("E", dim, dim).as_mutable()
    for i in range(dim):
        for j in range(dim):
            E[i,j] = sympy.Function(f"E_{i}_{j}")(*dofs)

    # Create an auxiliary symbol for the PK2 stress as a function depending on E(u,tetha)
    S_func = sympy.MatrixSymbol("S", dim, dim).as_mutable()
    for i in range(dim):
        for j in range(dim):
            S_func[i,j] = sympy.Function(f"S_{i}_{j}")(*(E.values()))

    # Substitute the previous auxiliary definition of the stress by S(E(u,tetha))
    SubstituteMatrixValue(rhs, S, S_func)

    # Set the substitution list to be applied after the differentiation
    n_dofs = rhs.shape[0]
    substitution_list = {}
    for i in range(dim):
        for j in range(dim):
            # First, for the dS/dE differentiation from the constitutive tensor symbols
            for k in range(dim):
                for l in range(dim):
                    t = sympy.diff(S_func[i, j], E[k, l])
                    substitution_list[t] = C[i, j, k, l]

            # Secondly, for the dE/dDOFs differentiation (note that we used the enriched strain)
            for i_dof in range(n_dofs):
                t = sympy.diff(E[i,j], dofs[i_dof])
                substitution_list[t] = sympy.diff(Emod_gauss[i,j], dofs[i_dof])

            # Finally, set the remaining E and S symbols as the enriched strain and the CL stress
            substitution_list[E[i,j]] = Emod_gauss[i,j]
            substitution_list[S_func[i,j]] = S[i,j]

    # Calculate the LHS from the previous derivatives
    n_dofs = rhs.shape[0]
    lhs = sympy.zeros(n_dofs, n_dofs)
    for i in range(n_dofs):
        for j in range(n_dofs):
                lhs[i,j] -= sympy.diff(rhs[i], dofs[j])
                lhs[i,j] = lhs[i,j].subs(substitution_list)
                if do_simplifications:
                    lhs[i, j] = sympy.simplify(lhs[i, j])
    lhs_out = OutputMatrix_CollectingFactors(w_g*lhs, "rLeftHandSideMatrix", mode, indentation_level=2, assignment_op="+=")

    # Calculate the output equivalent deformation gradient
    Fmod_gauss_out = OutputMatrix_CollectingFactors(Fmod_gauss, "r_eq_def_gradient", mode, indentation_level=1)
    det_Fmod_gauss_out = OutputScalar_CollectingFactors(Fmod_gauss.det(), "r_det_eq_def_gradient", mode, indentation_level=0)

    # Calculate the output equivalent strain
    Emod_gauss_out = OutputVector_CollectingFactors(StrainToVoigt(Emod_gauss), "r_eq_green_strain", mode, indentation_level=1)

    # Temporary save the automatic differentiation outputs
    tokens_dict = {
        f"substitute_rhs_{dim}D_{n_nodes}N" : rhs_out,
        f"substitute_lhs_{dim}D_{n_nodes}N" : lhs_out,
        f"substitute_def_gradient_{dim}D_{n_nodes}N" : Fmod_gauss_out,
        f"substitute_det_def_gradient_{dim}D_{n_nodes}N" : det_Fmod_gauss_out,
        f"substitute_green_strain_{dim}D_{n_nodes}N" : Emod_gauss_out}
    for token_key, token_val in tokens_dict.items():
        tokens_filenames.append(token_key)
        with open(f"{token_key}.sym",'w+') as token_file:
            token_file.write(token_val)

# Write the modified template
with open(template_filename, 'r') as template_file:
    with open(output_filename, 'w') as output_file:
        for line in template_file:
            has_token = False
            for token_key in tokens_filenames:
                if token_key in line:
                    has_token = True
                    with open(f"{token_key}.sym", 'r') as token_file:
                        output_file.writelines(token_file.readlines())
            if not has_token:
                output_file.write(line)

for token_key in tokens_filenames:
    kratos_utilities.DeleteFileIfExisting(f"{token_key}.sym")
