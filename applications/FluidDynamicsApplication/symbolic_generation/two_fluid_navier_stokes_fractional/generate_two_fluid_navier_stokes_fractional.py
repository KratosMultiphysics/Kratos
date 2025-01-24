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
dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
linearisation = "Picard"            # Iteration type. Options: "Picard", "FullNR"
divide_by_rho = True                # Divide by density in mass conservation equation
ASGS_stabilization = True           # Consider ASGS stabilization terms
mode = "c"                          # Output mode to a c++ file
time_integration="bdf2"             # TODO: Right now the formulation is validated only for BDF2 but there is a working version for generalized alpha.
adding_acceleration=True            # Add acceleration

if time_integration == "bdf2":
    output_filename = "two_fluid_navier_stokes_fractional.cpp"
    template_filename = "two_fluid_navier_stokes_fractional_template.cpp"
else:
    # TODO: There is a version with the generalised alpha time integration scheme but validation is requiered.
    err_msg = "Wrong time_integration. Given \'" + time_integration + "\'. Available options are \'bdf2\'."
    raise Exception(err_msg)

if (dim_to_compute == "2D"):
    dim_vector = [2]
elif (dim_to_compute == "3D"):
    dim_vector = [3]
elif (dim_to_compute == "Both"):
    dim_vector = [2,3]

## Read the template file
templatefile = open(template_filename)
outstring = templatefile.read()

for dim in dim_vector:

    if(dim == 2):
        nnodes = 3
        strain_size = 3
    else:
        nnodes = 4
        strain_size = 6

    impose_partion_of_unity = False
    gauss_weight = sympy.Symbol('w_gauss', positive = True)
    N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

    #define enrichment shape functions
    DNenr = DefineMatrix('DNenr',nnodes,dim)
    Nenr = DefineVector('Nenr',nnodes)

    ## Unknown fields definition
    v = DefineMatrix('v',nnodes,dim)            # Current step velocity (v(i,j) refers to velocity of node i component j)
    vn = DefineMatrix('vn',nnodes,dim)          # Previous step velocity
    vnn = DefineMatrix('vnn',nnodes,dim)        # 2 previous step velocity
    vnnn = DefineMatrix('vnnn', nnodes, dim)    # 3 previous step velocity TODO: We use this in case that the an is going to be obtained doing a BDF2 approximation
    p = DefineVector('p',nnodes)                # Pressure
    penr= DefineVector('penr',nnodes)	        # Enriched Pressure

    # Variables needed for the navier stokes fractional splitting
    vfrac = DefineMatrix('vfrac', nnodes, dim)  # Fractional velocity,

    ## Test functions definition
    w = DefineMatrix('w',nnodes,dim)            # Velocity field test function
    q = DefineVector('q',nnodes)                # Pressure field test function
    qenr = DefineVector('qenr' ,nnodes)	        # Enriched Pressure field test function

    ## Other data definitions
    f = DefineMatrix('f',nnodes,dim)            # Forcing term
    fn = DefineMatrix('fn',nnodes,dim)          # Previous step forcing term
    vmeshn = DefineMatrix('vmeshn',nnodes,dim)  # Previous step mesh velocity

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
    art_dyn_visc_coeff = sympy.Symbol('art_dyn_visc_coeff')

    ## Convective velocity definition
    if (linearisation == "Picard"):
        vconv = DefineMatrix('vconv',nnodes,dim)    # Convective velocity defined a symbol
    elif (linearisation == "FullNR"):
        vmesh = DefineMatrix('vmesh',nnodes,dim)    # Mesh velocity
        vconv = v - vmesh                           # Convective velocity defined as a velocity dependent variable

    vconv_gauss = vconv.transpose()*N              #Convective velocity on the current time step linearized
    vconv_old_gauss = vn.transpose()*N             #Convective velocity on the previous time step for the convective term of the fractional acceleration
    vconv_fractional = vfrac.transpose()*N         #Convective velocity corresponding to the convective fractional term.

    ## Compute the stabilization parameters
    vconv_gauss_norm = 0.0
    for i in range(0, dim):
        vconv_gauss_norm += vconv_gauss[i]**2
    vconv_gauss_norm = sympy.sqrt(vconv_gauss_norm)

    ## Data interpolation to the Gauss points
    ## Backward differences coefficients
    bdf0 = sympy.Symbol('bdf0')
    ## Part of the NS acceleration due to the fractional splitting
    acceleration = (bdf0*v -bdf0*vfrac)
    v_gauss = v.transpose()*N
    f_gauss = f.transpose()*N

    tau1 = 1.0/((rho*dyn_tau)/dt + (stab_c2*rho*vconv_gauss_norm)/h + (stab_c1*mu)/(h*h))   # Stabilization parameter 1
    tau2 = mu + (stab_c2*rho*vconv_gauss_norm*h)/stab_c1                                              # Stabilization parameter 2

    p_gauss = p.transpose()*N #NOTE: We evaluate p-related terms at n+1 as temporal component makes no sense in this case for both time integration schemes
    penr_gauss = penr.transpose()*Nenr
    w_gauss = w.transpose()*N
    q_gauss = q.transpose()*N
    qenr_gauss = qenr.transpose()*Nenr
    accel_gauss = acceleration.transpose()*N

    ## Gradients computation
    grad_v = DN.transpose()*v
    grad_v_fractional = DN.transpose()*vfrac
    grad_v_old = DN.transpose()*vn
    grad_w = DN.transpose()*w
    grad_q = DN.transpose()*q
    grad_qenr = DNenr.transpose()*qenr
    grad_p = DN.transpose()*p
    grad_penr = DNenr.transpose()*penr
    grad_sym_v_voigt = grad_sym_voigtform(DN, v)
    # Symmetric gradient of w in Voigt notation
    grad_sym_w_voigt = grad_sym_voigtform(DN, w)
    # Recall that the grad(w):sigma contraction equals grad_sym(w)*sigma in Voigt notation since sigma is a symmetric tensor.

    ## Divergence computation
    div_v = div(DN,v)
    div_w = div(DN,w)
    div_vconv = div(DN,vconv)

    # All the convective term definition
    convective_term = (vconv_gauss.transpose()*grad_v)                             #Convective term on the current time step NS equation
    convective_n_term = (vconv_old_gauss.transpose()*grad_v_old)                   #Convective term for the previous time step acceleration term (fractional splitting)
    convective_frac_term = vconv_fractional.transpose()*grad_v_fractional          #Fractional convective term

    # accel_n = bdf0*vn+bdf1*vnn+bdf2*vnnn #TODO: Approximating the acceleration with BDF2
    accel_n = (vn-vnn)/dt
    accel_gauss_n = accel_n.transpose()*N

    ## Galerkin Functional
    rv_galerkin = rho*w_gauss.transpose()*f_gauss - rho*w_gauss.transpose()*accel_gauss -rho*w_gauss.transpose()*convective_term.transpose() - grad_sym_w_voigt.transpose()*stress + div_w*p_gauss + rho*w_gauss.transpose()*convective_frac_term.transpose()

    if adding_acceleration:
        # Adding fractional acceleration convective part
        rv_galerkin -= rho*w_gauss.transpose()*convective_n_term.transpose()
        # Adding fractional acceleration partial part
        rv_galerkin -= rho*w_gauss.transpose()*accel_gauss_n

    if (divide_by_rho):
        rv_galerkin += q_gauss*(volume_error_ratio - div_v[0,0])
    else:
        rv_galerkin += rho*q_gauss*(volume_error_ratio - div_v[0,0])

    # Stabilization functional terms
    # Momentum conservation residual
    # Note that the viscous stress term is dropped since linear elements are used

    vel_residual = rho*f_gauss - rho*accel_gauss - rho*convective_term.transpose() -grad_p +rho*convective_frac_term.transpose()
    if adding_acceleration:
        vel_residual -=rho*(accel_gauss_n+convective_n_term.transpose())

    # Mass conservation residual
    if (divide_by_rho):
        mas_residual = -div_v[0,0] + volume_error_ratio
    else:
         mas_residual = -rho*div_v[0,0] + rho*volume_error_ratio

    vel_subscale = tau1*vel_residual
    mas_subscale = tau2*mas_residual

    # Compute the ASGS stabilization terms using the momentum and mass conservation residuals above
    if (divide_by_rho):
        rv_stab = grad_q.transpose()*vel_subscale
    else:
        rv_stab = rho*grad_q.transpose()*vel_subscale

    rv_stab += rho*vconv_gauss.transpose()*grad_w*vel_subscale
    rv_stab += rho*div_vconv*w_gauss.transpose()*vel_subscale
    rv_stab += div_w*mas_subscale

    ## Add the stabilization terms to the original residual terms
    if (ASGS_stabilization):
        rv = rv_galerkin + rv_stab
    else:
        rv = rv_galerkin

    ## Define DOFs and test function vectors
    dofs = sympy.zeros(nnodes*(dim+1), 1)
    testfunc = sympy.zeros(nnodes*(dim+1), 1)

    for i in range(0,nnodes):
        for k in range(0,dim):
            dofs[i*(dim+1)+k] = v[i,k]
            testfunc[i*(dim+1)+k] = w[i,k]
        dofs[i*(dim+1)+dim] = p[i,0]
        testfunc[i*(dim+1)+dim] = q[i,0]

    ## Compute LHS and RHS
    # For the RHS computation one wants the residual of the previous iteration (residual based formulation). By this reason the stress is
    # included as a symbolic variable, which is assumed to be passed as an argument from the previous iteration database.
    rhs = Compute_RHS(rv.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(gauss_weight*rhs, "rRHS", mode, assignment_op='+=')

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    # Note that the 'stress' (symbolic variable) is substituted by 'C*grad_sym_v_voigt' for the LHS differenctiation. Otherwise the velocity terms
    # within the velocity symmetryc gradient would not be considered in the differenctiation, meaning that the stress would be considered as
    # a velocity independent constant in the LHS.
    SubstituteMatrixValue(rhs, stress, C*grad_sym_v_voigt)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = OutputMatrix_CollectingFactors(gauss_weight * lhs, "rLHS", mode, assignment_op='+=')
    #Enrichment Functional
    ##  K V   x    =  b + rhs_eV
    ##  H Kee penr =  rhs_ee

    vel_residual_enr = rho*f_gauss - rho*(accel_gauss + convective_term.transpose()) - grad_p - grad_penr+rho*convective_frac_term.transpose()

    if adding_acceleration:
        vel_residual_enr -= rho*(accel_gauss_n+convective_n_term.transpose())

    vel_subscale_enr = vel_residual_enr * tau1
    rv_galerkin_enriched = div_w*penr_gauss

    if (divide_by_rho):
        rv_galerkin_enriched += qenr_gauss*(volume_error_ratio - div_v[0,0])
        rv_stab_enriched = grad_qenr.transpose()*vel_subscale_enr
        rv_stab_enriched -= grad_q.transpose()*tau1*grad_penr
    else:
        rv_galerkin_enriched += qenr_gauss*rho*(volume_error_ratio - div_v[0,0])
        rv_stab_enriched = rho*grad_qenr.transpose()*vel_subscale_enr
        rv_stab_enriched -= rho*grad_q.transpose()*tau1*grad_penr

    rv_stab_enriched -= rho*vconv_gauss.transpose()*grad_w*tau1*grad_penr
    rv_stab_enriched -= rho*div_vconv*w_gauss.transpose()*tau1*grad_penr
    rv_enriched = rv_galerkin_enriched
    ## Add the stabilization terms to the original residual terms
    if (ASGS_stabilization):
        rv_enriched += rv_stab_enriched

    dofs_enr=sympy.zeros(nnodes,1)
    testfunc_enr = sympy.zeros(nnodes,1)

    for i in range(0,nnodes):
        dofs_enr[i]=penr[i,0]
        testfunc_enr[i]=qenr[i,0]

    ##  K V   x    =  b + rhs_eV
    ##  H Kee penr =  rhs_ee
    rhs_eV, V   = Compute_RHS_and_LHS(rv_enriched, testfunc, dofs_enr, do_simplifications)
    rhs_ee, H   = Compute_RHS_and_LHS(rv_enriched, testfunc_enr, dofs, do_simplifications)
    rhs_ee, Kee = Compute_RHS_and_LHS(rv_enriched, testfunc_enr, dofs_enr, do_simplifications)

    V_out = OutputMatrix_CollectingFactors(gauss_weight*V,"rV",mode,assignment_op='+=')
    H_out = OutputMatrix_CollectingFactors(gauss_weight*H,"rH",mode,assignment_op='+=')
    Kee_out = OutputMatrix_CollectingFactors(gauss_weight*Kee,"rKee",mode,assignment_op='+=')
    rhs_ee_out = OutputVector_CollectingFactors(gauss_weight*rhs_ee,"rRHS_ee",mode,assignment_op='+=')

    #  Calculate artificial dynamic viscosity in each Gauss point
    vel_residual_norm = vel_residual.norm()
    grad_v_norm = grad_v.norm()
    artificial_mu = 0.5*h*art_dyn_visc_coeff*(vel_residual_norm/grad_v_norm)

    grad_v_norm_out = OutputScalar(grad_v_norm, "grad_v_norm", mode)
    artificial_mu_out = OutputScalar(artificial_mu, "artificial_mu", mode)

    #####################################################################
    #####################################################################

    if(dim == 2):
        outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)

        outstring = outstring.replace("//substitute_enrichment_V_2D", V_out)
        outstring = outstring.replace("//substitute_enrichment_H_2D", H_out)
        outstring = outstring.replace("//substitute_enrichment_Kee_2D", Kee_out)
        outstring = outstring.replace("//substitute_enrichment_rhs_ee_2D", rhs_ee_out)

        outstring = outstring.replace("//substitute_artificial_mu_2D_3N", artificial_mu_out)
        outstring = outstring.replace("//substitute_artificial_mu_grad_v_norm_2D_3N", grad_v_norm_out)

    elif(dim == 3):
        outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

        outstring = outstring.replace("//substitute_enrichment_V_3D", V_out)
        outstring = outstring.replace("//substitute_enrichment_H_3D", H_out)
        outstring = outstring.replace("//substitute_enrichment_Kee_3D", Kee_out)
        outstring = outstring.replace("//substitute_enrichment_rhs_ee_3D", rhs_ee_out)

        outstring = outstring.replace("//substitute_artificial_mu_3D_4N", artificial_mu_out)
        outstring = outstring.replace("//substitute_artificial_mu_grad_v_norm_3D_4N", grad_v_norm_out)

#We write in the file
out = open(output_filename,'w')
out.write(outstring)
out.close()

