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


## Symbolic generation settings
do_simplifications = False
dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
linearisation = "Picard"            # Iteration type. Options: "Picard", "FullNR"
ASGS_stabilization = True           # Consider ASGS stabilization terms
mode = "c"                          # Output mode to a c++ file
adding_acceleration = True         # Whether to add acceleration


output_filename = "two_fluid_navier_stokes_fractional_convection_element.cpp"
template_filename = "two_fluid_navier_stokes_fractional_convection_template.cpp"

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
    vfrac = DefineMatrix('vfrac', nnodes, dim)
    vn = DefineMatrix('vn',nnodes,dim)          # Previous step velocity
    vnn = DefineMatrix('vnn',nnodes,dim)        # 2 previous step velocity
    vmeshn = DefineMatrix('vmeshn',nnodes,dim)  # Previous step mesh velocity
    # vnnn = DefineMatrix('vnnn', nnodes, dim)  # 3 previous step velocity. In case of using an BDF2 old acceleration
    vconv = DefineMatrix('vconv',nnodes,dim)    # Convective velocity defined a symbol and in the element template is it substituted for vfrac from the previous iteration.

    ## Test functions definition.
    w = DefineMatrix('w',nnodes,dim)            # Velocity field test function

    # Scalar symbols.
    dt  = sympy.Symbol('dt', positive = True)
    nu  = sympy.Symbol('nu', positive = True)
    mu  = sympy.Symbol('mu', positive = True)
    tau1 = sympy.Symbol('tau1', positive = True)
    tau2 = sympy.Symbol('tau2', positive = True)
    h = sympy.Symbol('h', positive = True)
    dyn_tau = sympy.Symbol('dyn_tau', positive = True)
    stab_c1 = sympy.Symbol('stab_c1', positive = True)
    stab_c2 = sympy.Symbol('stab_c2', positive = True)
    bdf0 = sympy.Symbol('bdf0')
    bdf1 = sympy.Symbol('bdf1')
    bdf2 = sympy.Symbol('bdf2')
    acceleration = (bdf0*vfrac + bdf1*vn+bdf2*vnn) # This acceleration is calculated in a differnt way of a bdf2 since it is considered in the fractional splitting.
    vconv_gauss = vconv.transpose()*N
    vconv_gauss_norm = 0.0
    for i in range(0, dim):
        vconv_gauss_norm += vconv_gauss[i]**2
    vconv_gauss_norm = sympy.sqrt(vconv_gauss_norm)
    tau1 = 1.0/((dyn_tau/dt) + (stab_c2*vconv_gauss_norm)/h)

    # Interpolate to the gauss points.
    accel_gauss = acceleration.transpose()*N
    vfrac_gauss = vfrac.transpose()*N
    vn_gauss = vn.transpose()*N
    w_gauss = w.transpose()*N

    ## Gradients computation
    grad_vfrac = DN.transpose()*vfrac
    grad_v_old = DN.transpose()*vn
    grad_w = DN.transpose()*w

    # Convective term definition
    convective_term = (vconv_gauss.transpose()*grad_vfrac)

    # Accelerated convection definitoin
    # "an = dv_n/dt +v_n*grad_v_old"
    if adding_acceleration:
        # dv_d/dt    TODO:It is posible to use a BDF2
        accel_n = (vn-vnn)/dt
        accel_gauss_n = accel_n.transpose()*N
        # Convective past term definition vn*gradv_n
        convective_n_term = (vn_gauss.transpose()*grad_v_old)

    ## Galerkin Functional
    rv_galerkin =-w_gauss.transpose()*accel_gauss -w_gauss.transpose()*convective_term.transpose()

    # Adding acceleration
    if adding_acceleration:
        rv_galerkin+= w_gauss.transpose()*accel_gauss_n+ w_gauss.transpose()*convective_n_term.transpose()

    # Stabilization functional terms
    # Momentum conservation residual
    # Note that the viscous stress term is dropped since linear elements are used
    vel_residual = -accel_gauss - convective_term.transpose()
    if adding_acceleration:
        vel_residual+=  accel_gauss_n + convective_n_term.transpose()
    vel_subscale = tau1*vel_residual

    # Compute the ASGS stabilization terms using the momentum and mass conservation residuals above
    rv_stab =vconv_gauss.transpose()*grad_w*vel_subscale
    div_vconv = div(DN,vconv)
    rv_stab += div_vconv*w_gauss.transpose()*vel_subscale

    ## Add the stabilization terms to the original residual terms
    if (ASGS_stabilization):
        rv = rv_galerkin + rv_stab
    else:
        rv = rv_galerkin

    ## Define DOFs and test function vectors
    dofs = sympy.zeros(nnodes*(dim), 1)
    testfunc = sympy.zeros(nnodes*(dim), 1)

    for i in range(0,nnodes):
        for k in range(0,dim):
            dofs[i*(dim)+k] = vfrac[i,k]
            testfunc[i*(dim)+k] = w[i,k]

    ## Compute LHS and RHS
    # For the RHS computation one wants the residual of the previous iteration (residual based formulation).
    rhs = Compute_RHS(rv.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(gauss_weight*rhs, "rRHS", mode, assignment_op='+=')

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications)
    lhs_out = OutputMatrix_CollectingFactors(gauss_weight * lhs, "rLHS", mode, assignment_op='+=')

    #####################################################################
    #####################################################################

    if(dim == 2):
        outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)


    elif(dim == 3):
        outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

#We write in the file
out = open(output_filename,'w')
out.write(outstring)
out.close()

