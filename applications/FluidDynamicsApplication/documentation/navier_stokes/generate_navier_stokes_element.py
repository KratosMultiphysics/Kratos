from sympy import *
from KratosMultiphysics import *
from sympy_fe_utilities import *

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
# to the usual Navier-Stokes equations that are intended to act as a soft artificial compressibility, which is controlled by the value of "c".
# In the derivation of the residual equations the space variations of the density are considered to be close to 0.

## Symbolic generation settings
do_simplifications = False
dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
linearisation = "FullNR"            # Linearisation type. Options: "Picard", "FullNR"
divide_by_rho = True                # Divide by density in mass conservation equation
artificial_compressibility = True   # Consider an artificial compressibility
mode = "c"                          # Output mode to a c++ file

if (dim_to_compute == "2D"):
    dim_vector = [2]
elif (dim_to_compute == "3D"):
    dim_vector = [3]
elif (dim_to_compute == "Both"):
    dim_vector = [2,3]

## Read the template file
templatefile = open("navier_stokes_cpp_template.cpp")
outstring = templatefile.read()

for dim in dim_vector:

    if(dim == 2):
        nnodes = 3
        strain_size = 3
    elif(dim == 3):
        nnodes = 4
        strain_size = 6

    impose_partion_of_unity = False
    N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

    ## Unknown fields definition
    v = DefineMatrix('v',nnodes,dim)            # Current step velocity (v(i,j) refers to velocity of node i component j)
    vn = DefineMatrix('vn',nnodes,dim)          # Previous step velocity
    vnn = DefineMatrix('vnn',nnodes,dim)        # 2 previous step velocity
    p = DefineVector('p',nnodes)                # Pressure
    pn = DefineVector('pn',nnodes)              # Previous step pressure
    pnn = DefineVector('pnn',nnodes)            # 2 previous step pressure

    ## Test functions definition
    w = DefineMatrix('w',nnodes,dim)            # Velocity field test function
    q = DefineVector('q',nnodes)                # Pressure field test function

    ## Other data definitions
    f = DefineMatrix('f',nnodes,dim)            # Forcing term

    ## Constitutive matrix definition
    C = DefineSymmetricMatrix('C',strain_size,strain_size)

    ## Stress vector definition
    stress = DefineVector('stress',strain_size)

    ## Other simbols definition
    c   = Symbol('c',positive= True)            # Wave length number
    dt  = Symbol('dt', positive = True)         # Time increment
    rho = Symbol('rho', positive = True)        # Density
    nu  = Symbol('nu', positive = True)         # Kinematic viscosity (mu/rho)
    mu  = Symbol('mu', positive = True)         # Dynamic viscosity
    tau1 = Symbol('tau1', positive = True)      # Stabilization parameter 1
    tau2 = Symbol('tau2', positive = True)      # Stabilization parameter 2

    ## Backward differences coefficients
    bdf0 = Symbol('bdf0')
    bdf1 = Symbol('bdf1')
    bdf2 = Symbol('bdf2')

    ## Data interpolation to the Gauss points
    f_gauss = f.transpose()*N
    v_gauss = v.transpose()*N

    if (linearisation == "Picard"):
        vconv = DefineMatrix('vconv',nnodes,dim)    # Convective velocity defined a symbol
        vconv_gauss = vconv.transpose()*N
    elif (linearisation == "FullNR"):
        vmesh = DefineMatrix('vmesh',nnodes,dim)    # Mesh velocity
        vconv = v - vmesh                           # Convective velocity defined as a velocity dependent variable
        vconv_gauss = vconv.transpose()*N

    accel_gauss = (bdf0*v + bdf1*vn + bdf2*vnn).transpose()*N

    p_gauss = p.transpose()*N
    pder_gauss = (bdf0*p + bdf1*pn + bdf2*pnn).transpose()*N

    w_gauss = w.transpose()*N
    q_gauss = q.transpose()*N

    ## Gradients computation (fluid dynamics gradient)
    grad_w = DfjDxi(DN,w)
    grad_q = DfjDxi(DN,q)
    grad_p = DfjDxi(DN,p)
    grad_v = DfjDxi(DN,v)
    grad_vconv = DfjDxi(DN,vconv)

    div_w = div(DN,w)
    div_v = div(DN,v)
    div_vconv = div(DN,vconv)

    grad_sym_v = grad_sym_voigtform(DN,v)       # Symmetric gradient of v in Voigt notation
    grad_w_voigt = grad_sym_voigtform(DN,w)     # Symmetric gradient of w in Voigt notation -> TODO: CHECK THE *2 IN THE CROSS TERM
    # Recall that the grad(w):sigma contraction equals grad_sym(w)*sigma in Voigt notation since sigma is a symmetric tensor.

    # Convective term definition
    convective_term = (vconv_gauss.transpose()*grad_v)

    ## Compute galerkin functional
    # Navier-Stokes functional
    if (divide_by_rho):
        rv_galerkin = rho*w_gauss.transpose()*f_gauss - rho*w_gauss.transpose()*accel_gauss - rho*w_gauss.transpose()*convective_term.transpose() - grad_w_voigt.transpose()*stress + div_w*p_gauss - q_gauss*div_v
        if (artificial_compressibility):
            rv_galerkin -= (1/(rho*c*c))*q_gauss*pder_gauss
    else:
        rv_galerkin = rho*w_gauss.transpose()*f_gauss - rho*w_gauss.transpose()*accel_gauss - rho*w_gauss.transpose()*convective_term.transpose() - grad_w_voigt.transpose()*stress + div_w*p_gauss - rho*q_gauss*div_v
        if (artificial_compressibility):
            rv_galerkin -= (1/(c*c))*q_gauss*pder_gauss

    ##  Stabilization functional terms
    # Momentum conservation residual
    # Note that the viscous stress term is dropped since linear elements are used
    vel_residual = rho*f_gauss - rho*(accel_gauss + convective_term.transpose()) - grad_p

    # Mass conservation residual
    if (divide_by_rho):
        mas_residual = -div_v
        if (artificial_compressibility):
            mas_residual -= (1/(rho*c*c))*pder_gauss
    else:
        mas_residual = -rho*div_v
        if (artificial_compressibility):
            mas_residual -= (1/(c*c))*pder_gauss

    # Compute the ASGS stabilization terms using the momentum and mass conservation residuals above
    if (divide_by_rho):
        rv_stab = tau1*grad_q.transpose()*vel_residual
    else:
        rv_stab = tau1*rho*grad_q.transpose()*vel_residual
    rv_stab += tau1*rho*(vconv_gauss.transpose()*grad_w + div_vconv*w_gauss.transpose())*vel_residual
    rv_stab += tau2*div_w*mas_residual

    ## Add the stabilization terms to the original residual terms
    rv = rv_galerkin + rv_stab

    ## Define DOFs and test function vectors
    dofs = Matrix( zeros(nnodes*(dim+1), 1) )
    testfunc = Matrix( zeros(nnodes*(dim+1), 1) )

    for i in range(0,nnodes):

        # Velocity DOFs and test functions
        for k in range(0,dim):
            dofs[i*(dim+1)+k] = v[i,k]
            testfunc[i*(dim+1)+k] = w[i,k]

        # Pressure DOFs and test functions
        dofs[i*(dim+1)+dim] = p[i,0]
        testfunc[i*(dim+1)+dim] = q[i,0]

    ## Compute LHS and RHS
    # For the RHS computation one wants the residual of the previous iteration (residual based formulation). By this reason the stress is
    # included as a symbolic variable, which is assumed to be passed as an argument from the previous iteration database.
    rhs = Compute_RHS(rv.copy(), testfunc, False)
    rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    # Note that the 'stress' (symbolic variable) is substituted by 'C*grad_sym_v' for the LHS differenctiation. Otherwise the velocity terms
    # within the velocity symmetryc gradient would not be considered in the differenctiation, meaning that the stress would be considered as
    # a velocity independent constant in the LHS.
    SubstituteMatrixValue(rhs, stress, C*grad_sym_v)
    lhs = Compute_LHS(rhs, testfunc, dofs, False) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)


    if(dim == 2):
        outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)
    elif(dim == 3):
        outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

    ## Compute velocity subscale Gauss point value
    v_s_gauss = tau1*rho*(f_gauss - accel_gauss - convective_term.transpose()) - tau1*grad_p
    v_s_gauss_out = OutputVector_CollectingFactors(v_s_gauss, "v_s_gauss", mode)

    if(dim == 2):
        outstring = outstring.replace("//substitute_gausspt_subscale_2D", v_s_gauss_out)
    elif(dim == 3):
        outstring = outstring.replace("//substitute_gausspt_subscale_3D", v_s_gauss_out)

## Write the modified template
out = open("navier_stokes.cpp",'w')
out.write(outstring)
out.close()
