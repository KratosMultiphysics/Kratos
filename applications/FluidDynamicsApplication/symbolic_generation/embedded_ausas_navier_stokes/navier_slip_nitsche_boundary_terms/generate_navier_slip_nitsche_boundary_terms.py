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

## Symbolic generation settings
do_simplifications = False
dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
linearisation = "Picard"            # Iteration type. Options: "Picard", "FullNR"
mode = "c"                          # Output mode to a c++ file

## Set the dimensions to compute
if (dim_to_compute == "2D"):
    dim_vector = [2]
elif (dim_to_compute == "3D"):
    dim_vector = [3]
elif (dim_to_compute == "Both"):
    dim_vector = [2,3]

## Read the template file
templatefile = open("embedded_fluid_element_discontinuous_template.cpp")
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
    
    normal = DefineVector('normal',dim)         # Wall outwards normal vector
    g = DefineVector('g',dim)                   # Wall velocity (EMBEDDED_VELOCITY)

    ## Constitutive matrix definition
    C = DefineSymmetricMatrix('C',strain_size,strain_size)

    ## Stress vector definition
    stress = DefineVector('stress',strain_size)

    ## Other simbols definition
    c   = Symbol('c',positive= True)              # Wave length number
    dt  = Symbol('dt', positive = True)           # Time increment
    rho = Symbol('rho', positive = True)          # Density
    nu  = Symbol('nu', positive = True)           # Kinematic viscosity (mu/rho)
    mu  = Symbol('mu', positive = True)           # Dynamic viscosity
    h = Symbol('h', positive = True)              # Element size
    gamma = Symbol('gamma', positive = True)      # Penalty constant
    eps = Symbol('eps', positive = True)          # Epsilon (slip length)
    phi_u = Symbol('phi_u', positive = True)      # Nitsche stabilization constant
    dyn_tau = Symbol('dyn_tau', positive = True)  # Stabilization dynamic tau
    stab_c1 = Symbol('stab_c1', positive = True)  # Stabilization constant 1 
    stab_c2 = Symbol('stab_c2', positive = True)  # Stabilization constant 2
    adjoint = Symbol('adjoint')                   # Stabilization constant 2

    ## Backward differences coefficients
    bdf0 = Symbol('bdf0')
    bdf1 = Symbol('bdf1')
    bdf2 = Symbol('bdf2')

    ## Data interpolation to the Gauss points
    f_gauss = f.transpose()*N
    v_gauss = v.transpose()*N

    ## Convective velocity definition
    if (linearisation == "Picard"):
        vconv = DefineMatrix('vconv',nnodes,dim)    # Convective velocity defined a symbol
    elif (linearisation == "FullNR"):
        vmesh = DefineMatrix('vmesh',nnodes,dim)    # Mesh velocity
        vconv = v - vmesh                           # Convective velocity defined as a velocity dependent variable

    vconv_gauss = vconv.transpose()*N

    ## Compute the stabilization parameters
    vconv_gauss_norm = 0.0
    for i in range(0, dim):
        vconv_gauss_norm += vconv_gauss[i]**2
    vconv_gauss_norm = sqrt(vconv_gauss_norm)

    phi_u = mu + rho*vconv_gauss_norm*h + rho*h*h/dt                                        # Boundary term stabilization parameter

    ## Compute the rest of magnitudes at the Gauss points
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
    grad_sym_w = grad_sym_voigtform(DN,w)     # Symmetric gradient of w in Voigt notation
    # Recall that the grad(w):sigma contraction equals grad_sym(w)*sigma in Voigt notation since sigma is a symmetric tensor.

    # Define helper functions
    def VoigtStressToTensor(x):
        if x.shape[0] == 3:
            A = Matrix([ \
                [x[0], x[2]], \
                [x[2], x[1]]  \
                ])
        elif x.shape[0] == 6:
            A = Matrix([ \
                [x[0], x[3], x[5]], \
                [x[3], x[1], x[4]], \
                [x[5], x[4], x[2]]  \
                ])
        else:
            raise Exception("not supported")
        return A

    def VoigtStrainToTensor(x):
        if x.shape[0] == 3:
            A = Matrix([ \
                [x[0], x[2]/2], \
                [x[2]/2, x[1]]  \
                ])
        elif x.shape[0] == 6:
            A = Matrix([ \
                [x[0], x[3]/2, x[5]/2], \
                [x[3]/2, x[1], x[4]/2], \
                [x[5]/2, x[4]/2, x[2]]  \
                ])
        else:
            raise Exception("not supported")
        return A
        
    ## Define projection operators
    Pn = normal * normal.transpose()
    Pt = -Pn.copy()
    for i in range(dim):
        Pt[i,i] += 1

    ## Define nitsche terms
    ugPn = Pn*(v_gauss-g)
    ugPt = Pt*(v_gauss-g)

    nCgrad_sym_W = VoigtStressToTensor(C*grad_sym_w)*normal

    norm_imposition = ugPn
    tang_imposition = eps*Pt*VoigtStressToTensor(C*grad_sym_v)*normal + mu*ugPt

    tang_imposition_1 = eps * (VoigtStressToTensor(C*grad_sym_v)*normal).transpose() * Pt
    tang_imposition_2 = mu * (v_gauss - g).transpose() * Pt

    rv = \
        - (mu/(gamma*h)) * w_gauss.transpose() * norm_imposition \
        - (phi_u/(gamma*h)) * w_gauss.transpose() * norm_imposition \
        + norm_imposition.transpose() * normal * q_gauss \
        + adjoint * norm_imposition.transpose() * nCgrad_sym_W \
        - (1/(eps+gamma*h)) * tang_imposition_1 * w_gauss \
        - (1/(eps+gamma*h)) * tang_imposition_2 * w_gauss \
        + adjoint * ((gamma*h)/(eps+gamma*h)) * tang_imposition_1 * (2.0 * VoigtStrainToTensor(grad_sym_w) * normal) \
        + adjoint * ((gamma*h)/(eps+gamma*h)) * tang_imposition_2 * (2.0 * VoigtStrainToTensor(grad_sym_w) * normal) \
        # - adjoint * ((gamma*h)/(eps+gamma*h)) * tang_imposition_1 * (normal * q_gauss) \ # These term of the flux yields null contribution
        # - adjoint * ((gamma*h)/(eps+gamma*h)) * tang_imposition_2 * (normal * q_gauss) \ # These term of the flux yields null contribution

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
    rhs = Compute_RHS(rv.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    # Note that the 'stress' (symbolic variable) is substituted by 'C*grad_sym_v' for the LHS differenctiation. Otherwise the velocity terms
    # within the velocity symmetryc gradient would not be considered in the differenctiation, meaning that the stress would be considered as
    # a velocity independent constant in the LHS.
    SubstituteMatrixValue(rhs, stress, C*grad_sym_v)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)

    if(dim == 2):
        outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)
    elif(dim == 3):
        outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

## Write the modified template
out = open("embedded_fluid_element_discontinuous.cpp",'w')
out.write(outstring)
out.close()
