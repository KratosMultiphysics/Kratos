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
dim_to_compute = "2D"               # Spatial dimensions to compute. Options:  "2D","3D","Both"
linearisation = "Picard"            # Iteration type. Options: "Picard", "FullNR"
divide_by_rho = False               # Divide by density in mass conservation equation
artificial_compressibility = False  # Consider an artificial compressibility
FIC_stabilization = False           # Consider FIC stabilization terms. Not implemented
mode = "c"                          # Output mode to a c++ file

if (dim_to_compute == "2D"):
    dim_vector = [2]
    nvar = 1
elif (dim_to_compute == "3D"):
    dim_vector = [3]
elif (dim_to_compute == "Both"):
    dim_vector = [2,3]

## Read the template file
templatefile = open("helmholtz_cpp_template.cpp")
outstring = templatefile.read()

for dim in dim_vector:

    if(dim == 2):
        nnodes = 3
    elif(dim == 3):
        nnodes = 4
        strain_size = 6

    impose_partion_of_unity = False
    N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)  

    ## Unknown fields definition
    #~ v = DefineMatrix('v',nnodes,dim)            # Current step velocity (v(i,j) refers to velocity of node i component j)
    #~ vn = DefineMatrix('vn',nnodes,dim)          # Previous step velocity
    #~ vnn = DefineMatrix('vnn',nnodes,dim)        # 2 previous step velocity
    eta = DefineVector('eta',nnodes)            # Free surface elevation
    etan = DefineVector('etan',nnodes)          # Previous step free surface elevation
    etann = DefineVector('etann',nnodes)        # 2 previous step free surface elevation

    ## Test functions definition
    #~ w = DefineMatrix('w',nnodes,dim)            # Velocity field test function
    u = DefineVector('u',nnodes)                # Free surface elevation field test function

    ## Other data definitions
    #~ f = DefineMatrix('f',nnodes,dim)            # Forcing term

    ## Constitutive matrix definition
    #~ C = DefineSymmetricMatrix('C',strain_size,strain_size)

    ## Stress vector definition
    #~ stress = DefineVector('stress',strain_size)

    ## Other simbols definition
    dt = Symbol('dt', positive = True)          # Time increment
    H  = Symbol('H', positive = False)          # Bathymetry (H)
    h  = Symbol('h', positive = False)          # Depth (h)
    g  = Symbol('g', positive = True)           # gravity
    l  = Symbol('l', positive = True)           # Element size

    ## Backward differences coefficients
    bdf0 = Symbol('bdf0')
    bdf1 = Symbol('bdf1')
    bdf2 = Symbol('bdf2')

    ## Data interpolation to the Gauss points
    #~ f_gauss   = f.transpose()*N
    #~ v_gauss   = v.transpose()*N
    eta_gauss = eta.transpose()*N

    #~ ## Convective velocity definition
    #~ if (linearisation == "Picard"):
        #~ vconv = DefineMatrix('vconv',nnodes,dim)    # Convective velocity defined a symbol
    #~ elif (linearisation == "FullNR"):
        #~ vmesh = DefineMatrix('vmesh',nnodes,dim)    # Mesh velocity
        #~ vconv = v - vmesh                           # Convective velocity defined as a velocity dependent variable
#~ 
    #~ vconv_gauss = vconv.transpose()*N

 
    ## Compute the rest of magnitudes at the Gauss points
    #~ accel_gauss = (bdf0*v + bdf1*vn + bdf2*vnn).transpose()*N

    etadder_gauss = (bdf0*eta + bdf1*etan + bdf2*etann).transpose()*N

    #~ w_gauss = w.transpose()*N
    u_gauss = u.transpose()*N

    ## Gradients computation (fluid dynamics gradient)
    #~ grad_w = DfjDxi(DN,w)
    grad_u = DfjDxi(DN,u)
    grad_eta = DfjDxi(DN,eta)
    #~ grad_v = DfjDxi(DN,v)
    #~ grad_vconv = DfjDxi(DN,vconv)

    #~ div_w = div(DN,w)
    #~ div_v = div(DN,v)
    #~ div_vconv = div(DN,vconv)

    #~ grad_sym_v = grad_sym_voigtform(DN,v)       # Symmetric gradient of v in Voigt notation
    #~ grad_w_voigt = grad_sym_voigtform(DN,w)     # Symmetric gradient of w in Voigt notation
    #~ # Recall that the grad(w):sigma contraction equals grad_sym(w)*sigma in Voigt notation since sigma is a symmetric tensor.

    # Convective term definition
    #~ convective_term = (vconv_gauss.transpose()*grad_v)

    ## Compute galerkin functional
    # Helmholtz functional
    if (divide_by_rho):
        rv_galerkin = u_gauss.transpose()*etadder_gauss + grad_u.transpose()*g*H*grad_eta
        #~ rv_galerkin = rho*w_gauss.transpose()*f_gauss - rho*w_gauss.transpose()*accel_gauss - rho*w_gauss.transpose()*convective_term.transpose() - grad_w_voigt.transpose()*stress + div_w*p_gauss - q_gauss*div_v

    else:
        rv_galerkin = u_gauss.transpose()*etadder_gauss + grad_u.transpose()*g*H*grad_eta
        #~ rv_galerkin = rho*w_gauss.transpose()*f_gauss - rho*w_gauss.transpose()*accel_gauss - rho*w_gauss.transpose()*convective_term.transpose() - grad_w_voigt.transpose()*stress + div_w*p_gauss - rho*q_gauss*div_v


    #~ ##  Stabilization functional terms
    #~ # Momentum conservation residual
    #~ # Note that the viscous stress term is dropped since linear elements are used
    #~ vel_residual = rho*f_gauss - rho*accel_gauss -rho*convective_term.transpose() - grad_p
#~ 
    #~ # Mass conservation residual
    #~ if (divide_by_rho):
        #~ mas_residual = -div_v
        #~ if (artificial_compressibility):
            #~ mas_residual -= (1/(rho*c*c))*pder_gauss
    #~ else:
        #~ mas_residual = -rho*div_v
        #~ if (artificial_compressibility):
            #~ mas_residual -= (1/(c*c))*pder_gauss
#~ 
    #~ vel_subscale = tau1*vel_residual
    #~ mas_subscale = tau2*mas_residual
#~ 
    #~ # Compute the ASGS stabilization terms using the momentum and mass conservation residuals above
    #~ if (divide_by_rho):
        #~ rv_stab = grad_q.transpose()*vel_subscale
    #~ else:
        #~ rv_stab = rho*grad_q.transpose()*vel_subscale
    #~ rv_stab += rho*vconv_gauss.transpose()*grad_w*vel_subscale
    #~ rv_stab += rho*div_vconv*w_gauss.transpose()*vel_subscale
    #~ rv_stab += div_w*mas_subscale
#~ 
    #~ ## Add the stabilization terms to the original residual terms
    #~ if (ASGS_stabilization):
        #~ rv = rv_galerkin + rv_stab
    #~ else:
        #~ rv = rv_galerkin

    rv = rv_galerkin

    ## Define DOFs and test function vectors
    dofs = Matrix( zeros(nnodes*nvar, 1) )
    testfunc = Matrix( zeros(nnodes*nvar, 1) )

    for i in range(0,nnodes):

        #~ # Velocity DOFs and test functions
        #~ for k in range(0,dim):
            #~ dofs[i*nvar+k] = v[i,k]
            #~ testfunc[i*nvar+k] = w[i,k]

        # Free Surface elevation DOFs and test functions
        dofs[i*nvar] = eta[i,0]
        testfunc[i*nvar] = u[i,0]

    ## Compute LHS and RHS
    # For the RHS computation one wants the residual of the previous iteration (residual based formulation). By this reason the stress is
    # included as a symbolic variable, which is assumed to be passed as an argument from the previous iteration database.
    rhs = Compute_RHS(rv.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    # Note that the 'stress' (symbolic variable) is substituted by 'C*grad_sym_v' for the LHS differenctiation. Otherwise the velocity terms
    # within the velocity symmetryc gradient would not be considered in the differenctiation, meaning that the stress would be considered as
    # a velocity independent constant in the LHS.
    # SubstituteMatrixValue(rhs, stress, C*grad_sym_v)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)


    if(dim == 2):
        outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)
    elif(dim == 3):
        outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

    ## Compute velocity subscale Gauss point value
    #~ v_s_gauss = tau1*rho*(f_gauss - accel_gauss - convective_term.transpose()) - tau1*grad_p
    #~ v_s_gauss_out = OutputVector_CollectingFactors(v_s_gauss, "v_s_gauss", mode)
#~ 
    #~ if(dim == 2):
        #~ outstring = outstring.replace("//substitute_gausspt_subscale_2D", v_s_gauss_out)
    #~ elif(dim == 3):
        #~ outstring = outstring.replace("//substitute_gausspt_subscale_3D", v_s_gauss_out)

## Write the modified template
out = open("helmholtz.cpp",'w')
out.write(outstring)
out.close()
