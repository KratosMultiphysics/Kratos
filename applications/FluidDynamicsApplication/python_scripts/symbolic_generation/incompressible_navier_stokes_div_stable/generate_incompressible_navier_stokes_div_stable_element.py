import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Settings explanation
# DIMENSION TO COMPUTE:
# This symbolic generator is valid for both 2D and 3D cases. Since the element has been programed with a dimension template in Kratos,
# it is advised to set the dim_to_compute flag as "Both". In this case the generated .cpp file will contain both 2D and 3D implementations.
# LINEARISATION SETTINGS:
# FullNR considers the convective velocity as "v-vmesh", hence v is taken into account in the derivation of the LHS and RHS.
# Picard (a.k.a. QuasiNR) considers the convective velocity as "a", thus it is considered as a constant in the derivation of the LHS and RHS.
# DIVIDE BY RHO:
# If set to true, divides the mass conservation equation by rho in order to have a better conditioned matrix. Otherwise the original form is kept.
# ARTIFICIAL COMPRESSIBILITY:
# If set to true, the time derivative of the density is introduced in the mass conservation equation together with the state equation
# dp/drho=c^2 (being c the sound velocity). Besides, the velocity divergence is not considered to be 0. These assumptions add some extra terms
# to the usual Navier-Stokes equations that act as a weak compressibility controlled by the value of "c".
# CONVECTIVE TERM:
# If set to true, the convective term is taken into account in the calculation of the variational form. This allows generating both
# Navier-Stokes and Stokes elements.

## Symbolic generation settings
mode = "c"
do_simplifications = False
# dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
dim_to_compute = "2D"             # Spatial dimensions to compute. Options:  "2D","3D","Both"

convective_term = True
linearisation = "Picard" # Convective term linearisation type. Options: "Picard", "FullNR"
output_filename = "incompressible_navier_stokes_div_stable.cpp"
template_filename = "incompressible_navier_stokes_div_stable_cpp_template.cpp"

info_msg = "\n"
info_msg += "Element generator settings:\n"
info_msg += "\t - Dimension: " + dim_to_compute + "\n"
print(info_msg)

if dim_to_compute == "2D":
    dim_vector = [2]
    v_nodes_vector = [6] # tria
    p_nodes_vector = [3] # tria
elif dim_to_compute == "3D":
    dim_vector = [3]
    v_nodes_vector = [10] # tet
    p_nodes_vector = [4] # tet
elif dim_to_compute == "Both":
    dim_vector = [2, 3]
    v_nodes_vector = [6, 10] # tria, tet
    p_nodes_vector = [3, 4] # tria, tet

## Initialize the outstring to be filled with the template .cpp file
print("Reading template file \'"+ template_filename + "\'\n")
templatefile = open(template_filename)
outstring = templatefile.read()

for dim, v_n_nodes, p_n_nodes in zip(dim_vector, v_nodes_vector, p_nodes_vector):

    if dim == 2:
        strain_size = 3
    elif dim == 3:
        strain_size = 6

    N_v, DN_v = DefineShapeFunctions(v_n_nodes, dim, impose_partion_of_unity=False, shape_functions_name='N_v', gradients_name='DN_v')
    N_p, DN_p = DefineShapeFunctions(p_n_nodes, dim, impose_partion_of_unity=False, shape_functions_name='N_p', gradients_name='DN_p')

    ## Unknown fields definition
    v = DefineMatrix('r_v', v_n_nodes, dim)            # Current step velocity (v(i,j) refers to velocity of node i component j)
    vn = DefineMatrix('r_vn', v_n_nodes, dim)          # Previous step velocity
    vnn = DefineMatrix('r_vnn', v_n_nodes, dim)        # 2 previous step velocity
    p = DefineVector('r_p', p_n_nodes)                 # Pressure

    ## Fluid properties
    mu = sympy.Symbol('mu', positive = True)         # Dynamic viscosity
    rho = sympy.Symbol('rho', positive = True)       # Density

    ## Test functions definition
    w = DefineMatrix('w', v_n_nodes, dim)            # Velocity field test function
    q = DefineVector('q', p_n_nodes)                 # Pressure field test function

    ## Other data definitions
    f = DefineMatrix('r_f',v_n_nodes,dim)                 # Forcing term

    ## Constitutive matrix definition
    C = DefineSymmetricMatrix('C',strain_size,strain_size)

    ## Stress vector definition
    stress = DefineVector('r_stress',strain_size)

    ## Other simbols definition
    h = sympy.Symbol('h', positive = True)                        # Element characteristic size
    stab_c1 = sympy.Symbol('stab_c1', positive = True)            # Stabilization first constant
    stab_c2 = sympy.Symbol('stab_c2', positive = True)            # Stabilization second constant
    dyn_tau = sympy.Symbol('dyn_tau', positive = True)            # Stabilization dynamic tau
    dt = sympy.Symbol('rData.DeltaTime', positive = True)         # Time increment
    gauss_weight = sympy.Symbol('gauss_weight', positive = True)  # Integration point weight

    ## Backward differences coefficients
    bdf0 = sympy.Symbol('rData.BDF0')
    bdf1 = sympy.Symbol('rData.BDF1')
    bdf2 = sympy.Symbol('rData.BDF2')

    ## Convective velocity definition
    if convective_term:
        if (linearisation == "Picard"):
            vconv = DefineMatrix('vconv',v_n_nodes,dim)     # Convective velocity defined a symbol
        elif (linearisation == "FullNR"):
            vmesh = DefineMatrix('r_vmesh',v_n_nodes,dim)   # Mesh velocity
            vconv = v - vmesh                               # Convective velocity defined as a velocity dependent variable
        else:
            raise Exception("Wrong linearisation \'" + linearisation + "\' selected. Available options are \'Picard\' and \'FullNR\'.")
        vconv_gauss = vconv.transpose()*N_v

    ## Compute the rest of magnitudes at the Gauss points
    accel_gauss = (bdf0*v + bdf1*vn + bdf2*vnn).transpose()*N_v

    ## Data interpolation to the Gauss points
    f_gauss = f.transpose()*N_v

    v_gauss = v.transpose()*N_v
    p_gauss = p.transpose()*N_p

    w_gauss = w.transpose()*N_v
    q_gauss = q.transpose()*N_p

    ## Gradients computation (fluid dynamics gradient)
    grad_w = DfjDxi(DN_v, w)
    grad_q = DfjDxi(DN_p, q)

    grad_v = DfjDxi(DN_v,v)
    grad_p = DfjDxi(DN_p, p)

    div_w = div(DN_v,w)
    div_v = div(DN_v,v)

    if convective_term:
        div_vconv = div(DN_v, vconv)

    grad_sym_v = grad_sym_voigtform(DN_v,v)       # Symmetric gradient of v in Voigt notation
    grad_w_voigt = grad_sym_voigtform(DN_v,w)     # Symmetric gradient of w in Voigt notation
    # Recall that the grad(w):stress contraction equals grad_sym(w)*stress in Voigt notation since the stress is a symmetric tensor.

    # Convective term definition
    if convective_term:
        convective_term_gauss = (vconv_gauss.transpose()*grad_v)

    ## Compute galerkin functional
    # Navier-Stokes functional
    rv_galerkin = rho*w_gauss.transpose()*f_gauss - rho*w_gauss.transpose()*accel_gauss - grad_w_voigt.transpose()*stress + div_w*p_gauss - q_gauss*div_v
    if convective_term:
        rv_galerkin -= rho*w_gauss.transpose()*convective_term_gauss.transpose()

    ## Stabilization functional
    stab_norm_a = 0.0
    for i in range(dim):
        stab_norm_a += vconv_gauss[i]**2
    stab_norm_a = sympy.sqrt(stab_norm_a)
    tau = 1.0/(rho*dyn_tau/dt + stab_c2*rho*stab_norm_a/h + stab_c1*mu/h**2) # Stabilization operator

    mom_residual = rho*f_gauss - rho*accel_gauss -rho*convective_term_gauss.transpose() - grad_p #TODO: Here there should be the viscous component of the stress
    vel_subscale = tau * mom_residual

    rv_stab = grad_q.transpose()*vel_subscale
    rv_stab += rho*vconv_gauss.transpose()*grad_w*vel_subscale
    rv_stab += rho*div_vconv*w_gauss.transpose()*vel_subscale
    #TODO: Think about the viscous subscale component --> should we include it in the stress calculation?

    ## Define DOFs and test function vectors
    n_dofs = v_n_nodes * dim + p_n_nodes

    dofs = sympy.zeros(n_dofs, 1)
    testfunc = sympy.zeros(n_dofs, 1)

    # Velocity DOFs and test functions
    for i in range(v_n_nodes):
        for k in range(dim):
            dofs[i*dim + k] = v[i,k]
            testfunc[i*dim + k] = w[i,k]

    # Pressure DOFs and test functions
    for i in range(p_n_nodes):
        dofs[v_n_nodes*dim + i] = p[i,0]
        testfunc[v_n_nodes*dim + i] = q[i,0]

    ## Compute LHS and RHS
    # Add the stabilization to the Galerkin residual
    functional = rv_galerkin + rv_stab

    # For the RHS computation one wants the residual of the previous iteration (residual based formulation). By this reason the stress is
    # included as a symbolic variable, which is assumed to be passed as an argument from the previous iteration database.
    print(f"Computing {dim}D RHS Gauss point contribution\n")
    rhs = Compute_RHS(functional.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(gauss_weight*rhs, "rRHS", mode, assignment_op='+=')

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    # Note that the 'stress' (symbolic variable) is substituted by 'C*grad_sym_v' for the LHS differenctiation. Otherwise the velocity terms
    # within the velocity symmetryc gradient would not be considered in the differenctiation, meaning that the stress would be considered as
    # a velocity independent constant in the LHS.
    print(f"Computing {dim}D LHS Gauss point contribution\n")
    SubstituteMatrixValue(rhs, stress, C*grad_sym_v)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = OutputMatrix_CollectingFactors(gauss_weight*lhs, "rLHS", mode, assignment_op='+=')

    ## Replace the computed RHS and LHS in the template outstring
    outstring = outstring.replace(f"//substitute_lhs_{dim}D", lhs_out)
    outstring = outstring.replace(f"//substitute_rhs_{dim}D", rhs_out)

## Write the modified template
print("Writing output file \'" + output_filename + "\'")
out = open(output_filename,'w')
out.write(outstring)
out.close()
