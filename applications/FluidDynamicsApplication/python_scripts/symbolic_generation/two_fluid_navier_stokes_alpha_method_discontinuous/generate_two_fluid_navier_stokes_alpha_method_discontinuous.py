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


## Symbolic generation settings
do_simplifications = False
# dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
dim_to_compute = "2D"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
linearisation = "Picard"            # Iteration type. Options: "Picard", "FullNR"
ASGS_stabilization = True           # Consider ASGS stabilization terms
mode = "c"                          # Output mode to a c++ file
time_integration="alpha_method"

if time_integration == "alpha_method":
    output_filename = "two_fluid_navier_stokes_alpha_method_discontinuous.cpp"
    template_filename = "two_fluid_navier_stokes_alpha_method_discontinuous_template.cpp"
else:
    err_msg = f"Wrong time_integration. Given \'{time_integration}\'. Available option is \'alpha_method\'."
    raise Exception(err_msg)

if dim_to_compute == "2D":
    dim_vector = [2]
elif dim_to_compute == "3D":
    dim_vector = [3]
elif dim_to_compute == "Both":
    dim_vector = [2,3]

## Read the template file
templatefile = open(template_filename)
outstring = templatefile.read()

for dim in dim_vector:

    if(dim == 2):
        n_nodes = 3
        strain_size = 3
    else:
        n_nodes = 4
        strain_size = 6

    impose_partion_of_unity = False
    N, DN = DefineShapeFunctions(n_nodes, dim, impose_partion_of_unity)
    N_disc, DN_disc = DefineShapeFunctions(n_nodes, dim, impose_partion_of_unity)

    ## Define shape functions for p interpolation (standard)
    N = DefineVector('N', n_nodes)
    DN = DefineMatrix('DN', n_nodes, dim)

    ## Define shape functions for pressure enrichment
    N_enr_p = DefineVector('N_enr_p', n_nodes)
    DN_enr_p = DefineMatrix('DN_enr_p', n_nodes, dim)

    ## Define shape functions for velocity enrichment (bubble)
    # N_enr_vel = sympy.Symbol('N_enr_vel', positive=True)
    # DN_enr_vel = DefineVector('DN_enr_vel', dim)

    N_enr_vel_pos = sympy.Symbol('N_enr_vel_pos', positive=True)
    DN_enr_vel_pos = DefineVector('DN_enr_vel_pos', dim)

    N_enr_vel_neg = sympy.Symbol('N_enr_vel_neg', positive=True)
    DN_enr_vel_neg = DefineVector('DN_enr_vel_neg', dim)

    ## Unknown fields definition
    v = DefineMatrix('v', n_nodes, dim)         # Current step velocity (v(i,j) refers to velocity of node i component j)
    vn = DefineMatrix('vn', n_nodes, dim)       # Previous step velocity
    vnn = DefineMatrix('vnn', n_nodes, dim)     # 2 previous step velocity
    # v_enr = DefineVector('v_enr', dim)          # Velocity enrichment DOF
    v_enr_pos = DefineVector('v_enr_pos', dim)  # Velocity enrichment DOF
    v_enr_neg = DefineVector('v_enr_neg', dim)  # Velocity enrichment DOF
    p = DefineVector('p', n_nodes)              # Pressure
    p_enr= DefineVector('penr', n_nodes)	    # Enriched Pressure

    ## Test functions definition
    q = DefineVector('q', n_nodes)              # Pressure field test function
    w = DefineMatrix('w', n_nodes, dim)         # Velocity field test function
    q_enr = DefineVector('qenr' ,n_nodes)	    # Enriched pressure field test function
    # w_enr = DefineVector('w_enr', dim)          # Enriched velocity field test function
    w_enr_pos = DefineVector('w_enr_pos', dim)  # Enriched velocity field test function
    w_enr_neg = DefineVector('w_enr_neg', dim)  # Enriched velocity field test function

    ## Other data definitions
    f = DefineMatrix('f',n_nodes,dim)           # Forcing term
    fn = DefineMatrix('fn',n_nodes,dim)         # Previous step forcing term
    vmeshn = DefineMatrix('vmeshn',n_nodes,dim) # Previous step mesh velocity

    ## Constitutive matrix definition
    C = DefineSymmetricMatrix('C',strain_size,strain_size)

    ## Stress vector definition
    stress = DefineVector('stress',strain_size)

    ## Other simbols definition
    dt  = sympy.Symbol('dt', positive = True)
    rho = sympy.Symbol('rho', positive = True)
    nu  = sympy.Symbol('nu', positive = True)
    mu  = sympy.Symbol('mu', positive = True)
    tau1 = sympy.Symbol('tau1', positive = True)
    tau2 = sympy.Symbol('tau2', positive = True)
    h = sympy.Symbol('h', positive = True)
    dyn_tau = sympy.Symbol('dyn_tau', positive = True)
    stab_c1 = sympy.Symbol('stab_c1', positive = True)
    stab_c2 = sympy.Symbol('stab_c2', positive = True)
    volume_error_ratio = sympy.Symbol('volume_error_ratio')
    art_dyn_visc_coeff = sympy.Symbol('art_dyn_visc_coeff', positive = True)

    ## Gauss weight symbol definition
    # Note that it is required to output the LHS and RHS already multiplied by the Gauss to avoid the auxiliary assembly array
    gauss_weight = sympy.Symbol('gauss_weight', positive = True)

    ## Convective velocity definition
    if linearisation == "Picard":
        vconv = DefineMatrix('vconv', n_nodes, dim)    # Convective velocity defined a symbol
    elif linearisation == "FullNR":
        vmesh = DefineMatrix('vmesh', n_nodes, dim)    # Mesh velocity
        vconv = v - vmesh                              # Convective velocity defined as a velocity dependent variable
    vconv_gauss = vconv.transpose()*N

    ## Compute the stabilization parameters
    vconv_gauss_norm = 0.0
    for i in range(dim):
        vconv_gauss_norm += vconv_gauss[i]**2
    vconv_gauss_norm = sympy.sqrt(vconv_gauss_norm)

    ## Alpha method parameters
    max_sprectral_radius = sympy.Symbol('max_spectral_radius', positive = True)
    acceleration_alpha_method = DefineMatrix('acceleration_alpha_method', n_nodes, dim)
    alpha_m = 0.5*((3.0-max_sprectral_radius)/(1.0+max_sprectral_radius))
    alpha_f = 1.0/(1.0+max_sprectral_radius)
    gamma = 0.5 + alpha_m - alpha_f

    ## Alpha method affected variables
    f_alpha = fn + alpha_f*(f-fn)
    v_alpha = vn + alpha_f*(v-vn)
    v_gauss = v_alpha.transpose()*N
    f_gauss = f_alpha.transpose()*N
    acceleration_n = (v-vn)/(gamma*dt)+acceleration_alpha_method*(gamma-1.0)/gamma
    acceleration = acceleration_alpha_method + alpha_m*(acceleration_n-acceleration_alpha_method)

    ## Stabilization parameters
    tau1 = 1.0/((rho*dyn_tau)/dt + (stab_c2*rho*vconv_gauss_norm)/h + (stab_c1*mu)/(h*h))
    tau2 = mu + (stab_c2*rho*vconv_gauss_norm*h)/stab_c1

    ## Data interpolation to the Gauss points
    p_gauss = p.transpose()*N #NOTE: We evaluate p-related terms at n+1 as temporal component makes no sense in this case for both time integration schemes
    q_gauss = q.transpose()*N
    p_enr_gauss = p_enr.transpose()*N_enr_p
    q_enr_gauss = q_enr.transpose()*N_enr_p
    w_gauss = w.transpose()*N
    # w_enr_gauss = w_enr.transpose()*N_enr_vel
    w_enr_gauss = w_enr_pos.transpose()*N_enr_vel_pos + w_enr_neg.transpose()*N_enr_vel_neg
    accel_gauss = acceleration.transpose()*N

    ## Gradients computation
    grad_p = DN.transpose()*p
    grad_q = DN.transpose()*q

    grad_p_enr = DN_enr_p.transpose()*p_enr
    grad_q_enr = DN_enr_p.transpose()*q_enr

    grad_v = DN.transpose()*v_alpha
    grad_w = DN.transpose()*w

    # grad_v_enr = DN_enr_vel*v_enr.transpose()
    # grad_w_enr = DN_enr_vel*w_enr.transpose()
    grad_v_enr = DN_enr_vel_pos*v_enr_pos.transpose() + DN_enr_vel_neg*v_enr_neg.transpose()
    grad_w_enr = DN_enr_vel_pos*w_enr_pos.transpose() + DN_enr_vel_neg*w_enr_neg.transpose()

    div_v = div(DN,v_alpha)
    div_v_stabilization=div(DN,v)

    # div_v_enr = DN_enr_vel.transpose()*v_enr
    # div_w_enr = DN_enr_vel.transpose()*w_enr
    div_v_enr = DN_enr_vel_pos.transpose()*v_enr_pos + DN_enr_vel_neg.transpose()*v_enr_neg
    div_w_enr = DN_enr_vel_pos.transpose()*w_enr_pos + DN_enr_vel_neg.transpose()*w_enr_neg

    div_w = div(DN,w)
    div_vconv = div(DN,vconv)

    grad_sym_v_voigt = grad_sym_voigtform(DN, v_alpha) # Symmetric gradient of v in Voigt notation
    grad_sym_w_voigt = grad_sym_voigtform(DN, w)       # Symmetric gradient of w in Voigt notation

    grad_sym_v_enr_voigt = sympy.simplify(StrainToVoigt(0.5*(grad_v_enr+grad_v_enr.transpose()))) # Symmetric gradient of enriched v in Voigt notation
    grad_sym_w_enr_voigt = sympy.simplify(StrainToVoigt(0.5*(grad_w_enr+grad_w_enr.transpose()))) # Symmetric gradient of enriched w in Voigt notation

    # Convective term definition
    convective_term = vconv_gauss.transpose()*grad_v
    convective_term_enr = vconv_gauss.transpose()*grad_v_enr #TODO: Should we include v_enr in the convective velocity? Maybe save the one from prev. it and use it?

    ## Galerkin Functional
    rv_galerkin = rho*w_gauss.transpose()*f_gauss - rho*w_gauss.transpose()*accel_gauss - rho*w_gauss.transpose()*convective_term.transpose() - grad_sym_w_voigt.transpose()*stress + div_w*p_gauss
    rv_galerkin += q_gauss*(volume_error_ratio - div_v[0,0])

    ## Stabilization functional terms
    # Momentum conservation residual
    # Note that the viscous stress term is dropped since linear elements are used
    vel_residual = rho*f_gauss - rho*accel_gauss - rho*convective_term.transpose() - grad_p

    # Mass conservation residual
    mas_residual = volume_error_ratio - div_v_stabilization[0,0] #NOTE: The divergence of v in here is not computed at alpha

    ## Subscales calculation
    vel_subscale = tau1*vel_residual
    mas_subscale = tau2*mas_residual

    ## Compute the ASGS stabilization terms using the momentum and mass conservation residuals above
    rv_stab = grad_q.transpose()*vel_subscale
    rv_stab += rho*vconv_gauss.transpose()*grad_w*vel_subscale
    rv_stab += rho*div_vconv*w_gauss.transpose()*vel_subscale
    rv_stab += div_w*mas_subscale

    ## Add the stabilization terms to the original residual terms
    if ASGS_stabilization:
        rv = rv_galerkin + rv_stab
    else:
        rv = rv_galerkin

    ## Define DOFs and test function vectors
    dofs = sympy.zeros(n_nodes*(dim+1), 1)
    test_func = sympy.zeros(n_nodes*(dim+1), 1)

    for i in range(n_nodes):
        for k in range(dim):
            dofs[i*(dim+1)+k] = v[i,k]
            test_func[i*(dim+1)+k] = w[i,k]
        dofs[i*(dim+1)+dim] = p[i,0]
        test_func[i*(dim+1)+dim] = q[i,0]

    ## Compute LHS and RHS
    # For the RHS computation one wants the residual of the previous iteration (residual based formulation). By this reason the stress is
    # included as a symbolic variable, which is assumed to be passed as an argument from the previous iteration database.
    rhs = Compute_RHS(rv.copy(), test_func, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(gauss_weight*rhs, "rRHS", mode, indentation_level=2, assignment_op="+=")

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    # Note that the 'stress' (symbolic variable) is substituted by 'C*grad_sym_v_voigt' for the LHS differenctiation. Otherwise the velocity terms
    # within the velocity symmetryc gradient would not be considered in the differenctiation, meaning that the stress would be considered as
    # a velocity independent constant in the LHS.
    SubstituteMatrixValue(rhs, stress, C*grad_sym_v_voigt)
    lhs = Compute_LHS(rhs, test_func, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = OutputMatrix_CollectingFactors(gauss_weight*lhs, "rLHS", mode, indentation_level=2, assignment_op="+=")

    # Enrichment Functional
    ##  K V   x    =  b + rhs_eV
    ##  H Kee x_enr =  rhs_ee

    # Calculate subscales for the enrichment terms functional
    vel_residual_enr = - rho*convective_term_enr.transpose() - grad_p_enr
    mas_residual_enr = - div_v_enr[0,0] #NOTE: The divergence of v in here is not computed at alpha

    vel_subscale_enr = tau1 * vel_residual_enr
    mas_subscale_enr = tau2 * mas_residual_enr

    # Calculate the enrichment functional
    rv_galerkin_enriched = rho*w_enr_gauss*f_gauss - rho*w_enr_gauss*accel_gauss
    rv_galerkin_enriched += - rho*w_gauss.transpose()*convective_term_enr.transpose() - rho*w_enr_gauss*convective_term.transpose() - rho*w_enr_gauss*convective_term_enr.transpose()
    rv_galerkin_enriched += - grad_sym_w_enr_voigt.transpose()*stress
    rv_galerkin_enriched += + div_w_enr[0,0]*p_gauss + div_w*p_enr_gauss + div_w_enr[0,0]*p_enr_gauss
    rv_galerkin_enriched += q_enr_gauss*(volume_error_ratio - div_v[0,0] - div_v_enr[0,0]) - q_gauss*div_v_enr[0,0]

    # Calculate the stabilization functional
    rv_stab_enriched = (rho*div_vconv)*w_gauss.transpose()*vel_subscale_enr
    rv_stab_enriched += rho*vconv_gauss.transpose()*grad_w*vel_subscale_enr
    rv_stab_enriched += div_w*mas_subscale_enr
    rv_stab_enriched += grad_q.transpose()*vel_subscale_enr
    rv_stab_enriched += (rho*div_vconv)*w_enr_gauss*(vel_subscale + vel_subscale_enr)
    rv_stab_enriched += rho*vconv_gauss.transpose()*grad_w*(vel_subscale + vel_subscale_enr)
    rv_stab_enriched += div_w_enr*(mas_subscale + mas_subscale_enr)
    rv_stab_enriched += grad_q_enr.transpose()*(vel_subscale + vel_subscale_enr)

    # Add the stabilization terms to the original residual terms
    rv_enriched = rv_galerkin_enriched
    if ASGS_stabilization:
        rv_enriched += rv_stab_enriched

    # Set up the enrichment DOFs and test functions arrays
    # dofs_enr = sympy.zeros(n_nodes + dim,1)
    # test_func_enr = sympy.zeros(n_nodes + dim,1)
    dofs_enr = sympy.zeros(n_nodes + 2*dim, 1)
    test_func_enr = sympy.zeros(n_nodes + 2*dim, 1)

    for i in range(n_nodes):
        dofs_enr[i] = p_enr[i,0]
        test_func_enr[i] = q_enr[i,0]
    # for i in range(dim):
    #     dofs_enr[n_nodes + i] = v_enr[i,0]
    #     test_func_enr[n_nodes + i] = w_enr[i,0]
    for i in range(dim):
        dofs_enr[n_nodes + i] = v_enr_pos[i,0]
        test_func_enr[n_nodes + i] = w_enr_pos[i,0]
    for i in range(dim):
        dofs_enr[n_nodes + dim + i] = v_enr_neg[i,0]
        test_func_enr[n_nodes + dim + i] = w_enr_neg[i,0]

    ##  K V   x    =  b + rhs_eV
    ##  H Kee p_enr =  rhs_ee
    rhs_eV = Compute_RHS(rv_enriched.copy(), test_func, do_simplifications)
    SubstituteMatrixValue(rhs_eV, stress, C*grad_sym_v_enr_voigt) # Note that the non-enriched stress is already added to the standard RHS
    V = Compute_LHS(rhs_eV, test_func, dofs_enr, do_simplifications)

    rhs_ee = Compute_RHS(rv_enriched.copy(), test_func_enr, do_simplifications)
    SubstituteMatrixValue(rhs_eV, stress, C*(grad_sym_v_voigt + grad_sym_v_enr_voigt))
    H = Compute_LHS(rhs_ee, test_func_enr, dofs, do_simplifications)
    Kee = Compute_LHS(rhs_ee, test_func_enr, dofs_enr, do_simplifications)

    V_out = OutputMatrix_CollectingFactors(gauss_weight * V, "rV", mode, indentation_level=2, assignment_op="+=")
    H_out = OutputMatrix_CollectingFactors(gauss_weight * H, "rH", mode, indentation_level=2, assignment_op="+=")
    Kee_out = OutputMatrix_CollectingFactors(gauss_weight * Kee, "rKee", mode, indentation_level=2, assignment_op="+=")
    rhs_ee_out = OutputVector_CollectingFactors(gauss_weight * rhs_ee, "rRHSee", mode, indentation_level=2, assignment_op="+=")

    # ## Calculate artificial dynamic viscosity in each Gauss point
    # grad_v_norm = grad_v.norm()
    # vel_residual_norm = vel_residual.norm()
    # artificial_mu = 0.5*h*art_dyn_visc_coeff*(vel_residual_norm/grad_v_norm)

    # grad_v_norm_out = OutputScalar(grad_v_norm, "grad_v_norm", mode)
    # artificial_mu_out = OutputScalar(artificial_mu, "artificial_mu", mode)

    #####################################################################
    #####################################################################

    if dim == 2:
        outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)

        outstring = outstring.replace("//substitute_enrichment_V_2D", V_out)
        outstring = outstring.replace("//substitute_enrichment_H_2D", H_out)
        outstring = outstring.replace("//substitute_enrichment_Kee_2D", Kee_out)
        outstring = outstring.replace("//substitute_enrichment_rhs_ee_2D", rhs_ee_out)

        # outstring = outstring.replace("//substitute_artificial_mu_2D_3N", artificial_mu_out)
        # outstring = outstring.replace("//substitute_artificial_mu_grad_v_norm_2D_3N", grad_v_norm_out)

    elif dim == 3:
        outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

        outstring = outstring.replace("//substitute_enrichment_V_3D", V_out)
        outstring = outstring.replace("//substitute_enrichment_H_3D", H_out)
        outstring = outstring.replace("//substitute_enrichment_Kee_3D", Kee_out)
        outstring = outstring.replace("//substitute_enrichment_rhs_ee_3D", rhs_ee_out)

        # outstring = outstring.replace("//substitute_artificial_mu_3D_4N", artificial_mu_out)
        # outstring = outstring.replace("//substitute_artificial_mu_grad_v_norm_3D_4N", grad_v_norm_out)

#We write in the file
out = open(output_filename,'w')
out.write(outstring)
out.close()

