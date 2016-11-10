from sympy import *
from KratosMultiphysics import *
from sympy_fe_utilities import *

do_simplifications = False
dim_to_compute = "Both"       # Spatial dimensions to compute. Options:  "2D","3D","Both"
mode = "c"                    # Output mode to a c++ file

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
    vconv = DefineMatrix('vconv',nnodes,dim)    # Convective velocity
    p = DefineVector('p',nnodes)                # Pressure

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
    dt  = Symbol('dt', positive = True)
    rho = Symbol('rho', positive = True)
    nu  = Symbol('nu', positive = True)
    mu  = Symbol('mu', positive = True)
    tau1 = Symbol('tau1', positive = True)
    tau2 = Symbol('tau2', positive = True)

    ## Backward differences coefficients
    bdf0 = Symbol('bdf0')
    bdf1 = Symbol('bdf1')
    bdf2 = Symbol('bdf2')

    ## Data interpolation to the Gauss points
    f_gauss = f.transpose()*N
    v_gauss = v.transpose()*N
    vconv_gauss = vconv.transpose()*N
    accel_gauss = (bdf0*v+bdf1*vn+bdf2*vnn).transpose()*N
    p_gauss = p.transpose()*N

    w_gauss = w.transpose()*N 
    q_gauss = q.transpose()*N

    ## Gradients computation
    grad_v = DN.transpose()*v
    grad_w = DN.transpose()*w
    grad_q = DN.transpose()*q
    grad_p = DN.transpose()*p

    #~ B = MatrixB(DN)

    grad_sym_f = grad_sym_voigtform(DN,f)
    grad_sym_v = grad_sym_voigtform(DN,v)

    div_v = div(DN,v)
    div_w = div(DN,w)

    grad_w_voigt = grad_sym_voigtform(DN,w)

    # Convective term definition
    convective_term = (vconv_gauss.transpose()*grad_v)

    ## Compute galerkin functional
    # Navier-Stokes functional
    rv_galerkin = rho*w_gauss.transpose()*f_gauss - rho*w_gauss.transpose()*accel_gauss - rho*convective_term*w_gauss - grad_w_voigt.transpose()*stress + div_w*p_gauss - rho*q_gauss*div_v

    # Stabilization functional terms
    vel_residual = rho*f_gauss - rho*(accel_gauss + convective_term.transpose()) - grad_p

    rv_stab = tau1*rho*grad_q.transpose()*vel_residual
    rv_stab += tau1*rho*vel_residual.transpose()*grad_w*vconv_gauss
    rv_stab -= tau2*rho*div_w*div_v

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
    rhs,lhs = Compute_RHS_and_LHS(rv.copy(), testfunc, dofs, False)

    lhs_out = OutputMatrix_CollectingFactors(lhs,"lhs", mode)
    rhs_out = OutputVector_CollectingFactors(rhs,"rhs", mode) 

    if(dim == 2):
        outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)
    elif(dim == 3):
        outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

## Write the modified template
out = open("navier_stokes.cpp",'w')
out.write(outstring)
out.close()

