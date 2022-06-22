import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

# Symbolic generation settings
mode = "c"
dim_vect = [2]
# dim_vect = [2, 3]
do_simplifications = False
output_filename = "total_lagrangian_mixed_detJ_element.cpp"
template_filename = "total_lagrangian_mixed_detJ_element_template.cpp"
outstring = open(template_filename).read()

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
    invF_gauss = F_gauss.inv()
    cofF_gauss = j_gauss*(invF_gauss.transpose())

    # Calculate the strain tensors
    Fbar_gauss = (1/(j_gauss**(1/3)))*F_gauss # Deviatoric deformation gradient tensor #TODO: Decide about the dim
    # Fbar_gauss = (1/(j_gauss**(1/dim)))*F_gauss # Deviatoric deformation gradient tensor
    Cbar_gauss = Fbar_gauss.transpose() * Fbar_gauss # Deviatoric right Cauchy-Green strain tensor
    Ebar_gauss = 0.5*(Cbar_gauss - sympy.eye(dim,dim)) # Deviatoric Green strain tensor

    # Calculate the equivalent strain tensors
    Fmod_gauss = ((1.0 + th_gauss)**(1/3))*Fbar_gauss # Equivalent (enriched) deformation gradient tensor #TODO: Decide about the dim
    # Fmod_gauss = ((1.0 + th_gauss)**(1/dim))*Fbar_gauss # Equivalent (enriched) deformation gradient tensor
    Cmod_gauss = Fmod_gauss.transpose() * Fmod_gauss # Equivalent (enriched) right Cauchy-Green strain tensor
    Emod_gauss = 0.5*(Cmod_gauss - sympy.eye(dim,dim)) # Equivalent (enriched) Green strain tensor

    # Variational form
    tmp = (DoubleContraction(C, F_gauss.transpose()*F_gauss)).tomatrix()

    mom_first = DoubleContraction(grad_w_gauss, F_gauss* S)
    mom_second = (w_gauss.transpose() * b_gauss)[0]
    mom_aux_scalar = (tau_th / 3.0) * ((1+th_gauss)**(-1.0/3.0)) * (1.0/(j_gauss**(2.0/3.0))) #TODO: Decide about the dim
    # mom_aux_scalar = (tau_th / dim) * ((1+th_gauss)**((2.0-dim)/dim)) * (1.0/(j_gauss**(2/dim)))
    mom_stab = DoubleContraction(grad_w_gauss, mom_aux_scalar * (1 + th_gauss - j_gauss) * tmp)

    mass_first = q_gauss[0] * (1.0+th_gauss - j_gauss)
    mass_aux_scalar = (tau_u / 3.0) * ((j_gauss / (1.0+th_gauss))**((1.0)/3.0)) #TODO: Decide about the dim
    # mass_aux_scalar = (tau_u / dim) * ((j_gauss / (1.0+th_gauss))**((dim-2.0)/dim))
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
    rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)

    # Compute LHS (RHS(residual) differentiation w.r.t. the DOFs)
    # Note that 'S' (stress symbolic variable) is substituted by 'C*E' for the LHS differentiation, being E the enriched strain.
    # Otherwise the displacement terms inside the modified Green strain would not be considered in the differentiation
    SubstituteMatrixValue(rhs, S, DoubleContraction(C, Emod_gauss).tomatrix())
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications)
    lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)

    # Replace the computed RHS and LHS in the template outstring
    outstring = outstring.replace(f"//substitute_rhs_{dim}D_{n_nodes}N", rhs_out)
    outstring = outstring.replace(f"//substitute_lhs_{dim}D_{n_nodes}N", lhs_out)

    # Replace the equivalent deformation gradient in the template outstring
    Fmod_gauss_out = OutputMatrix_CollectingFactors(Fmod_gauss, "r_eq_def_gradient", mode)
    det_Fmod_gauss_out = OutputScalar_CollectingFactors(Fmod_gauss.det(), "r_det_eq_def_gradient", mode)
    outstring = outstring.replace(f"//substitute_def_gradient_{dim}D_{n_nodes}N", Fmod_gauss_out)
    outstring = outstring.replace(f"//substitute_det_def_gradient_{dim}D_{n_nodes}N", det_Fmod_gauss_out)

    # Replace the equivalent strain in the template outstring
    Emod_gauss_out = OutputVector_CollectingFactors(StrainToVoigt(Emod_gauss), "r_eq_green_strain", mode)
    outstring = outstring.replace(f"//substitute_green_strain_{dim}D_{n_nodes}N", Emod_gauss_out)

# Write the modified template
with open(output_filename, 'w') as f:
    f.write(outstring)
