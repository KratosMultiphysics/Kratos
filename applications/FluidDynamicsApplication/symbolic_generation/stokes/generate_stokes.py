import sympy
from KratosMultiphysics import sympy_fe_utilities as sfu


## Settings explanation
# DIMENSION TO COMPUTE:
# This symbolic generator is valid for both 2D and 3D cases. Since the element has been programed with a dimension template in Kratos,
# it is advised to set the dim_to_compute flag as "Both". In this case the generated .cpp file will contain both 2D and 3D implementations.
# DIVIDE BY RHO:
# If set to true, divides the mass conservation equation by rho in order to have a better conditioned matrix. Otherwise the original form is kept.

## Symbolic generation settings
do_simplifications = False
dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
divide_by_rho = True                # Divide by density in mass conservation equation
ASGS_stabilization = True           # Consider ASGS stabilization terms
mode = "c"                          # Output mode to a c++ file

# all linear elements
if (dim_to_compute == "2D"):
    dim_vector = [2, 2]
    nnodes_vector = [3, 4] # tria
elif (dim_to_compute == "3D"):
    dim_vector = [3, 3, 3]
    nnodes_vector = [4, 6, 8] # tet, prism, hex
elif (dim_to_compute == "Both"):
    dim_vector = [2, 2, 3, 3, 3]
    nnodes_vector = [3, 4, 4, 6, 8] # tria, quad, tet, prism, hex

## Read the template file
templatefile = open("symbolic_stokes_template.cpp")
outstring = templatefile.read()

for dim, nnodes in zip(dim_vector, nnodes_vector):
    if(dim == 2):
        strain_size = 3
    else:
        strain_size = 6

    impose_partion_of_unity = False
    N,DN = sfu.DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)
    ## Unknown fields definition
    v = sfu.DefineMatrix('v',nnodes,dim)            # Current step velocity (v(i,j) refers to velocity of node i component j)
    vn = sfu.DefineMatrix('vn',nnodes,dim)          # Previous step velocity
    vnn = sfu.DefineMatrix('vnn',nnodes,dim)        # 2 previous step velocity
    p = sfu.DefineVector('p',nnodes)                # Pressure

    ## Test functions definition
    w = sfu.DefineMatrix('w',nnodes,dim)            # Velocity field test function
    q = sfu.DefineVector('q',nnodes)                # Pressure field test function

    ## Other data definitions
    f = sfu.DefineMatrix('f',nnodes,dim)            # Forcing term

    ## Constitutive matrix definition
    C = sfu.DefineSymmetricMatrix('C',strain_size,strain_size)

    ## Stress vector definition
    stress = sfu.DefineVector('stress',strain_size)

    ## Other simbols definition
    dt  = sfu.Symbol('dt', positive = True)
    dt_inv  = sfu.Symbol('dt_inv', positive = True)
    rho = sfu.Symbol('rho', positive = True)
    nu  = sfu.Symbol('nu', positive = True)
    mu  = sfu.Symbol('mu', positive = True)
    tau1 = sfu.Symbol('tau1', positive = True)
    tau2 = sfu.Symbol('tau2', positive = True)
    h = sfu.Symbol('h', positive = True)
    dyn_tau = sfu.Symbol('dyn_tau', positive = True)
    stab_c1 = sfu.Symbol('stab_c1', positive = True)
    ##stab_c2 = sfu.Symbol('stab_c2', positive = True)

    ## Backward differences coefficients
    bdf0 = sfu.Symbol('bdf0')
    bdf1 = sfu.Symbol('bdf1')
    bdf2 = sfu.Symbol('bdf2')


    ## Data interpolation to the Gauss points
    f_gauss = f.transpose()*N
    v_gauss = v.transpose()*N

    tau1 = 1.0/((rho*dyn_tau)*dt_inv + (stab_c1*mu)/(h*h))   # Stabilization parameter 1
    tau2 = mu # Stabilization parameter 2


    ## Values at Gauss point
    accel_gauss = (bdf0*v + bdf1*vn + bdf2*vnn).transpose()*N
    p_gauss = p.transpose()*N
    w_gauss = w.transpose()*N
    q_gauss = q.transpose()*N

    ## Gradients computation
    grad_v = DN.transpose()*v
    grad_w = DN.transpose()*w
    grad_q = DN.transpose()*q
    grad_p = DN.transpose()*p

    div_v = sfu.div(DN,v)
    div_w = sfu.div(DN,w)


    grad_sym_v_voigt = sfu.grad_sym_voigtform(DN,v)     # Symmetric gradient of v in Voigt notation
    grad_sym_w_voigt = sfu.grad_sym_voigtform(DN,w)     # Symmetric gradient of w in Voigt notation

    ## Galerkin Functional
    rv_galerkin = rho*w_gauss.transpose()*f_gauss - rho*w_gauss.transpose()*accel_gauss - grad_sym_w_voigt.transpose()*stress + div_w*p_gauss

    if (divide_by_rho):
        rv_galerkin -= q_gauss*div_v
    else:
        rv_galerkin -= rho*q_gauss*div_v

    # Stabilization functional terms
    # Momentum conservation residual
    # Note that the viscous stress term is dropped since linear elements are used
    vel_residual = rho*f_gauss - rho*accel_gauss - grad_p

    # Mass conservation residual
    if (divide_by_rho):
        mas_residual = -div_v
    else:
        mas_residual = -rho*div_v

    vel_subscale = tau1*vel_residual
    mas_subscale = tau2*mas_residual

    # Compute the ASGS stabilization terms using the momentum and mass conservation residuals above
    if (divide_by_rho):
        rv_stab = grad_q.transpose()*vel_subscale
    else:
        rv_stab = rho*grad_q.transpose()*vel_subscale

    rv_stab += div_w*mas_subscale

    ## Add the stabilization terms to the original residual terms
    if (ASGS_stabilization):
        rv = rv_galerkin + rv_stab
    else:
        rv = rv_galerkin

    ## Define DOFs and test function vectors
    dofs = sfu.Matrix( sfu.zeros(nnodes*(dim+1), 1) )
    testfunc = sfu.Matrix( sfu.zeros(nnodes*(dim+1), 1) )

    for i in range(0,nnodes):
        for k in range(0,dim):
            dofs[i*(dim+1)+k] = v[i,k]
            testfunc[i*(dim+1)+k] = w[i,k]
        dofs[i*(dim+1)+dim] = p[i,0]
        testfunc[i*(dim+1)+dim] = q[i,0]


    ## Compute LHS and RHS
    # For the RHS computation one wants the residual of the previous iteration (residual based formulation). By this reason the stress is
    # included as a symbolic variable, which is assumed to be passed as an argument from the previous iteration database.
    rhs = sfu.Compute_RHS(rv.copy(), testfunc, do_simplifications)
    rhs_out = sfu.OutputVector_CollectingFactors(rhs, "rhs", mode)

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    # Note that the 'stress' (symbolic variable) is substituted by 'C*grad_sym_v_voigt' for the LHS differenctiation. Otherwise the velocity terms
    # within the velocity symmetryc gradient would not be considered in the differenctiation, meaning that the stress would be considered as
    # a velocity independent constant in the LHS.
    sfu.SubstituteMatrixValue(rhs, stress, C*grad_sym_v_voigt)
    lhs = sfu.Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = sfu.OutputMatrix_CollectingFactors(lhs, "lhs", mode)

    #####################################################################
    #####################################################################

    if(dim == 2 and nnodes == 3):
        outstring = outstring.replace("//substitute_lhs_2D3", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D3", rhs_out)
    elif(dim == 2 and nnodes == 4):
        outstring = outstring.replace("//substitute_lhs_2D4", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D4", rhs_out)
    elif(dim == 3 and nnodes == 4):
        outstring = outstring.replace("//substitute_lhs_3D4", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D4", rhs_out)
    elif(dim == 3 and nnodes == 6):
        outstring = outstring.replace("//substitute_lhs_3D6", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D6", rhs_out)
    elif(dim == 3 and nnodes == 8):
        outstring = outstring.replace("//substitute_lhs_3D8", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D8", rhs_out)

#We write in the file
out = open("symbolic_stokes.cpp",'w')
out.write(outstring)
out.close()

