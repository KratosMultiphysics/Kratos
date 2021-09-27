from sympy import *
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

##TO DO: include (div_vconv) term ?? -> commented now

## Symbolic generation settings
do_simplifications = False
dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
linearisation = "Picard"            # Iteration type. Options: "Picard", "FullNR"
divide_by_rho = True                # Divide by density in mass conservation equation
ASGS_stabilization = True           # Consider ASGS stabilization terms
mode = "c"                          # Output mode to a c++ file


if (dim_to_compute == "2D"):
    dim_vector = [2]
elif (dim_to_compute == "3D"):
    dim_vector = [3]
elif (dim_to_compute == "Both"):
    dim_vector = [2,3]

## Read the template file
templatefile = open("two_fluid_navier_stokes_template_CN.cpp")
outstring = templatefile.read()

for dim in dim_vector:

    if(dim == 2):
        nnodes = 3
        strain_size = 3
    else:
        nnodes = 4
        strain_size = 6

    impose_partion_of_unity = False
    N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

    #define enrichment shape functions
    DNenr = DefineMatrix('DNenr',nnodes,dim)
    Nenr = DefineVector('Nenr',nnodes)

    theta=Symbol('theta',positive=True)

    ## Unknown fields definition
    v = DefineMatrix('v',nnodes,dim)            # Current step velocity (v(i,j) refers to velocity of node i component j)
    vn = DefineMatrix('vn',nnodes,dim)          # Previous step velocity
    v_cn = theta*v+(1-theta)*vn
    p = DefineVector('p',nnodes)                # Pressure
    pn = DefineVector('pn',nnodes)              # Previous steo pressure
    p_cn = theta*p+(1-theta)*pn
    #FIXME: Enriched Pressure hacemos algo ??        # Enriched Pressure
    penr_cn= DefineVector('penr_cn',nnodes)
    ## Test functions definition
    w = DefineMatrix('w',nnodes,dim)            # Velocity field test function
    q = DefineVector('q',nnodes)                # Pressure field test function
    qenr = DefineVector('qenr' ,nnodes)	        # Enriched Pressure field test function

    ## Other data definitions
    f = DefineMatrix('f',nnodes,dim)            # Forcing term

    ## Constitutive matrix definition
    C = DefineSymmetricMatrix('C',strain_size,strain_size)

    ## Stress vector definition
    stress = DefineVector('stress',strain_size)

    ## Other simbols definition
    dt  = Symbol('dt', positive = True)
    rho = Symbol('rho', positive = True)
    nu  = Symbol('nu', positive = True)
    mu  = Symbol('mu', positive = True)
    tau1 = Symbol('tau1', positive = True)
    tau2 = Symbol('tau2', positive = True)
    h = Symbol('h', positive = True)
    dyn_tau = Symbol('dyn_tau', positive = True)

    stab_c1 = Symbol('stab_c1', positive = True)
    stab_c2 = Symbol('stab_c2', positive = True)
    K_darcy = Symbol('K_darcy', positive = True)
    volume_error_ratio = Symbol('volume_error_ratio')

    ## Data interpolation to the Gauss points
    f_gauss = f.transpose()*N
    v_cn_gauss=v_cn.transpose()*N

    ## Convective velocity definition
    ## FIXME: vconv 0.5 evaluated
    if (linearisation == "Picard"):
        vconv_cn = DefineMatrix('vconv_CN',nnodes,dim)    # Convective velocity defined a symbol
    elif (linearisation == "FullNR"):
        vmesh = DefineMatrix('vmesh',nnodes,dim)          # Mesh velocity
        vconv_cn = v_cn - vmesh                           # Convective velocity defined as a velocity dependent variable

    vconv_cn_gauss = vconv_cn.transpose()*N



    ## Compute the stabilization parameters
    vconv_cn_gauss_norm = 0.0
    for i in range(0, dim):
        vconv_cn_gauss_norm += vconv_cn_gauss[i]**2
    vconv_cn_gauss_norm = sqrt(vconv_cn_gauss_norm)


    tau1 = 1.0/((rho*dyn_tau)/(dt) + (stab_c2*rho*vconv_cn_gauss_norm)/h + (stab_c1*mu)/(h*h) + K_darcy)   # Stabilization parameter 1
    tau2 = mu + (stab_c2*rho*vconv_cn_gauss_norm*h)/stab_c1                                    # Stabilization parameter 2


    ## Values at Gauss point
    ## FIXME: LA ACELERACION DEBE SER OTRA ??????
    accel_gauss_cn=((v-vn)/dt).transpose()*N
    p_cn_gauss=p_cn.transpose()*N
    # FIXME: LA ENRICHMENTE ???
    penr_cn_gauss = penr_cn.transpose()*Nenr
    w_gauss = w.transpose()*N
    q_gauss = q.transpose()*N
    qenr_gauss = qenr.transpose()*Nenr

    ## Gradients computation
    #FIXME: GRADIENTE ???
    grad_v_cn=DN.transpose()*v_cn
    grad_w = DN.transpose()*w
    grad_q = DN.transpose()*q
    grad_qenr = DNenr.transpose()*qenr
    grad_p_cn=DN.transpose()*p_cn
    grad_penr_cn = DNenr.transpose()*penr_cn
    # FIXME: DIVERGENCIA??
    div_v_cn=div(DN,v_cn)
    div_w = div(DN,w)
    # div_vconv = div(DN,vconv_cn)
    # div_vconv_cn=div(DN,vconv)

    # FIXME: GRADIENTE??
    grad_sym_v_cn_voigt = grad_sym_voigtform(DN,v_cn)   # Symmetric gradient of v in Voigt notation
    grad_sym_w_voigt = grad_sym_voigtform(DN,w)     # Symmetric gradient of w in Voigt notation
    # Recall that the grad(w):sigma contraction equals grad_sym(w)*sigma in Voigt notation since sigma is a symmetric tensor.


    # Convective term definition
    convective_cn_term = (vconv_cn_gauss.transpose()*grad_v_cn)

    ## Galerkin Functional
    rv_galerkin_cn = rho*w_gauss.transpose()*f_gauss - rho*w_gauss.transpose()*accel_gauss_cn - rho*w_gauss.transpose()*convective_cn_term.transpose() - grad_sym_w_voigt.transpose()*stress + div_w*p_cn_gauss
    rv_galerkin_cn -= w_gauss.transpose()*K_darcy*v_cn_gauss #Darcy Term

    if (divide_by_rho):
        rv_galerkin_cn += q_gauss*(volume_error_ratio - div_v_cn[0,0])
    else:
        rv_galerkin_cn += rho*q_gauss*(volume_error_ratio - div_v_cn[0,0])

    # Stabilization functional terms
    # Momentum conservation residual
    # Note that the viscous stress term is dropped since linear elements are used

    vel_residual_cn = rho*f_gauss - rho*accel_gauss_cn - rho*convective_cn_term.transpose() - grad_p_cn - K_darcy*v_cn_gauss

    # Mass conservation residual
    if (divide_by_rho):
        mas_residual_cn = -div_v_cn[0,0] + volume_error_ratio
    else:
        mas_residual_cn = -rho*div_v_cn[0,0] + rho*volume_error_ratio

    vel_subscale_cn = tau1*vel_residual_cn
    mas_subscale_cn = tau2*mas_residual_cn

    # Compute the ASGS stabilization terms using the momentum and mass conservation residuals above

    if (divide_by_rho):
        rv_stab_cn = grad_q.transpose()*vel_subscale_cn
    else:
        rv_stab_cn = rho*grad_q.transpose()*vel_subscale_cn


    rv_stab_cn += rho*vconv_cn_gauss.transpose()*grad_w*vel_subscale_cn
    #rv_stab += rho*div_vconv*w_gauss.transpose()*vel_subscale
    rv_stab_cn += div_w*mas_subscale_cn
    rv_stab_cn -= w_gauss.transpose()*K_darcy*vel_subscale_cn


    if (ASGS_stabilization):
        rv_cn = rv_galerkin_cn + rv_stab_cn
    else:
        rv_cn = rv_galerkin_cn


    ## Define DOFs and test function vectors
    dofs = Matrix( zeros(nnodes*(dim+1), 1) )
    testfunc = Matrix( zeros(nnodes*(dim+1), 1) )

    for i in range(0,nnodes):
        for k in range(0,dim):
            dofs[i*(dim+1)+k] = v[i,k]
            testfunc[i*(dim+1)+k] = w[i,k]
        dofs[i*(dim+1)+dim] = p[i,0]
        testfunc[i*(dim+1)+dim] = q[i,0]


    ## Compute LHS and RHS
    # For the RHS computation one wants the residual of the previous iteration (residual based formulation). By this reason the stress is
    # included as a symbolic variable, which is assumed to be passed as an argument from the previous iteration database.
    rhs = Compute_RHS(rv_cn.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    # Note that the 'stress' (symbolic variable) is substituted by 'C*grad_sym_v_cn_voigt' for the LHS differenctiation. Otherwise the velocity terms
    # within the velocity symmetryc gradient would not be considered in the differenctiation, meaning that the stress would be considered as
    # a velocity independent constant in the LHS.
    SubstituteMatrixValue(rhs, stress, C*grad_sym_v_cn_voigt)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)

    #Enrichment Functional
    ##  K V   x    =  b + rhs_eV
    ##  H Kee penr =  rhs_ee

    vel_residual_enr_cn = rho*f_gauss - rho*(accel_gauss_cn + convective_cn_term.transpose()) - grad_p_cn - K_darcy*v_cn_gauss - grad_penr_cn
    vel_subscale_enr_cn = vel_residual_enr_cn * tau1

    #FIXME: Enrichmente_without cn
    rv_galerkin_enriched_cn = div_w*penr_cn_gauss
    # rv_galerkin_enriched = div_w*penr_gauss

    if (divide_by_rho):
        rv_galerkin_enriched_cn += qenr_gauss*(volume_error_ratio - div_v_cn[0,0])
        rv_stab_enriched_cn = grad_qenr.transpose()*vel_subscale_enr_cn
        rv_stab_enriched_cn -= grad_q.transpose()*tau1*grad_penr_cn
    else:
        rv_galerkin_enriched_cn += qenr_gauss*rho*(volume_error_ratio - div_v_cn[0,0])
        rv_stab_enriched_cn = rho*grad_qenr.transpose()*vel_subscale_enr_cn
        rv_stab_enriched_cn -= rho*grad_q.transpose()*tau1*grad_penr_cn

    rv_stab_enriched_cn -= rho*vconv_cn_gauss.transpose()*grad_w*tau1*grad_penr_cn
    rv_stab_enriched_cn += w_gauss.transpose()*K_darcy*tau1*grad_penr_cn
    #rv_stab_enriched -= rho*div_vconv*w_gauss.transpose()*tau1*grad_penr

    rv_enriched_cn = rv_galerkin_enriched_cn
    ## Add the stabilization terms to the original residual terms
    if (ASGS_stabilization):
        rv_enriched_cn += rv_stab_enriched_cn

    dofs_enr=Matrix( zeros(nnodes,1) )
    testfunc_enr = Matrix( zeros(nnodes,1) )

    for i in range(0,nnodes):
        dofs_enr[i]=penr_cn[i,0]
        testfunc_enr[i]=qenr[i,0]

    ##  K V   x    =  b + rhs_eV
    ##  H Kee penr =  rhs_ee
    rhs_eV, V   = Compute_RHS_and_LHS(rv_enriched_cn, testfunc, dofs_enr, do_simplifications)
    rhs_ee, H   = Compute_RHS_and_LHS(rv_enriched_cn, testfunc_enr, dofs, do_simplifications)
    rhs_ee, Kee = Compute_RHS_and_LHS(rv_enriched_cn, testfunc_enr, dofs_enr, do_simplifications)

    V_out = OutputMatrix_CollectingFactors(V,"V",mode)
    H_out = OutputMatrix_CollectingFactors(H,"H",mode)
    Kee_out = OutputMatrix_CollectingFactors(Kee,"Kee",mode)
    rhs_ee_out = OutputVector_CollectingFactors(rhs_ee,"rhs_ee",mode)






    #####################################################################
    #####################################################################

    if(dim == 2):
        outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)

        outstring = outstring.replace("//substitute_enrichment_V_2D", V_out)
        outstring = outstring.replace("//substitute_enrichment_H_2D", H_out)
        outstring = outstring.replace("//substitute_enrichment_Kee_2D", Kee_out)
        outstring = outstring.replace("//substitute_enrichment_rhs_ee_2D", rhs_ee_out)

    elif(dim == 3):
        outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

        outstring = outstring.replace("//substitute_enrichment_V_3D", V_out)
        outstring = outstring.replace("//substitute_enrichment_H_3D", H_out)
        outstring = outstring.replace("//substitute_enrichment_Kee_3D", Kee_out)
        outstring = outstring.replace("//substitute_enrichment_rhs_ee_3D", rhs_ee_out)

#We write in the file
out = open("two_fluid_navier_stokes_CN.cpp",'w')
out.write(outstring)
out.close()

