import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Symbolic generation settings
mode = "c"
divide_by_rho = True                # Divide the mass conservation equation by rho
ASGS_stabilization = True           # Consider ASGS stabilization terms
add_incompressibility_error = False     # Add the incompressibility error to the viscous stress response

do_simplifications = False
output_filename = "axisymmetric_navier_stokes.cpp"
template_filename = "axisymmetric_navier_stokes_cpp_template.cpp"

info_msg = f"\n"
info_msg += f"Element generator settings:\n"
info_msg += f"\t - ASGS stabilization: {ASGS_stabilization}\n"
info_msg += f"\t - Divide mass conservation by rho: {divide_by_rho}\n"
info_msg += f"\t - Add incompressibility error to viscous stress: {add_incompressibility_error}\n"
print(info_msg)

dim_vector = [2, 2]
n_nodes_vector = [3, 4] # tria, quad

## Initialize the outstring to be filled with the template .cpp file
print("Reading template file \'"+ template_filename + "\'\n")
templatefile = open(template_filename)
outstring = templatefile.read()

for dim, n_nodes in zip(dim_vector, n_nodes_vector):
    ## Define shape functions, shape function gradients and integration weight
    impose_partion_of_unity = False
    gauss_weight = sympy.Symbol('w_gauss', positive = True) # Integration weight defined at cpp
    N,DN = DefineShapeFunctions(n_nodes, dim, impose_partion_of_unity)

    ## Unknown fields definition
    v = DefineMatrix('v',n_nodes,dim)            # Current step velocity (v(i,j) refers to velocity of node i component j)
    v_n = DefineMatrix('v_n',n_nodes,dim)        # Previous step velocity
    v_nn = DefineMatrix('v_nn',n_nodes,dim)      # 2 previous step velocity
    p = DefineVector('p',n_nodes)                # Pressure

    ## Fluid properties
    mu = sympy.Symbol('mu', positive = True) # Dynamic viscosity
    rho = sympy.Symbol('rho', positive = True) # Density

    ## Test functions definition
    w = DefineMatrix('w',n_nodes,dim)            # Velocity field test function
    q = DefineVector('q',n_nodes)                # Pressure field test function

    ## Other data definitions
    f = DefineMatrix('f',n_nodes,dim)            # Forcing term

    ## Other symbol definitions
    y = sympy.Symbol('y', positive = True)      # y-coordinate (radius)
    h = sympy.Symbol('h', positive = True)      # Element size
    dt  = sympy.Symbol('dt', positive = True)   # Time increment
    dyn_tau = sympy.Symbol('dyn_tau', positive = True)
    stab_c1 = sympy.Symbol('stab_c1', positive = True)
    stab_c2 = sympy.Symbol('stab_c2', positive = True)

    ## Backward differences coefficients
    bdf0 = sympy.Symbol('bdf0')
    bdf1 = sympy.Symbol('bdf1')
    bdf2 = sympy.Symbol('bdf2')

    ## Data interpolation to the Gauss points
    f_gauss = f.transpose()*N
    v_gauss = v.transpose()*N
    p_gauss = p.transpose()*N
    w_gauss = w.transpose()*N
    q_gauss = q.transpose()*N

    ## Convective velocity definition
    v_conv = DefineMatrix('v_conv', n_nodes, dim) # Convective velocity defined a symbol (this is required to avoid its differentiation)
    v_conv_gauss = v_conv.transpose()*N

    ## Compute the stabilization parameters
    stab_norm_a = 0.0
    for i in range(0, dim):
        stab_norm_a += v_conv_gauss[i]**2
    stab_norm_a = sympy.sqrt(stab_norm_a)
    tau1 = 1.0/((rho*dyn_tau)/dt + (stab_c2*rho*stab_norm_a)/h + (stab_c1*mu)/(h*h)) # Stabilization parameter 1
    tau2 = mu + (stab_c2*rho*stab_norm_a*h)/stab_c1                                  # Stabilization parameter 2

    ## Compute the rest of magnitudes at the Gauss points
    accel_gauss = (bdf0*v + bdf1*v_n + bdf2*v_nn).transpose()*N

    ## Gradients computation (fluid dynamics gradient)
    grad_w = DfjDxi(DN,w)
    grad_q = DfjDxi(DN,q)
    grad_p = DfjDxi(DN,p)
    grad_v = DfjDxi(DN,v)
    grad_v_conv = DfjDxi(DN,v_conv)

    ## Compute galerkin functional
    # Navier-Stokes functional
    galerkin_functional = rho*w_gauss[0]*f_gauss[0] - rho*w_gauss[0]*accel_gauss[0] - rho*w_gauss[0]*v_conv_gauss[1]*grad_v[1,0] - rho*w_gauss[0]*v_conv_gauss[0]*grad_v[0,0] + grad_w[0,0]*p_gauss[0] - mu*grad_w[1,0]*grad_v[1,0] - mu*grad_w[0,0]*grad_v[0,0]
    galerkin_functional += rho*w_gauss[1]*f_gauss[1] - rho*w_gauss[1]*accel_gauss[1] - rho*w_gauss[1]*v_conv_gauss[1]*grad_v[1,1] - rho*w_gauss[1]*v_conv_gauss[0]*grad_v[0,1] + (w_gauss[1]/y + grad_w[1,1])*p_gauss[0] - mu*grad_w[1,1]*grad_v[1,1] - mu*grad_w[0,1]*grad_v[0,1]
    if add_incompressibility_error:
        galerkin_functional -= 2.0*mu*grad_w[0,0]*grad_v[0,0]
        galerkin_functional -= 2.0*mu*(grad_w[1,1]*grad_v[1,1] + w_gauss[1]*v_gauss[1]/y**2)
    if divide_by_rho:
        galerkin_functional += - q_gauss[0]*v_gauss[1]/y - q_gauss[0]*grad_v[1,1] - q_gauss[0]*grad_v[0,0]
    else:
        galerkin_functional += - rho*q_gauss[0]*v_gauss[1]/y - rho*q_gauss[0]*grad_v[1,1] - rho*q_gauss[0]*grad_v[0,0]

    ##  Stabilization functional terms
    # Momentum conservation residual
    # Note that the second order viscous stress terms are dropped since linear elements are used
    vel_residual_x = rho*f_gauss[0] - rho*accel_gauss[0] - rho*v_conv_gauss[1]*grad_v[1,0] - rho*v_conv_gauss[0]*grad_v[0,0] - grad_p[0] + (mu/y)*grad_v[1,0]
    vel_residual_y = rho*f_gauss[1] - rho*accel_gauss[1] - rho*v_conv_gauss[1]*grad_v[1,1] - rho*v_conv_gauss[0]*grad_v[0,1] - grad_p[1] - mu*v_gauss[1]/y**2 + (mu/y)*grad_v[1,1]

    # Mass conservation residual
    if divide_by_rho:
        mas_residual = - v_gauss[1]/y - grad_v[1,1] - grad_v[0,0]
    else:
        mas_residual = - rho*v_gauss[1]/y - rho*grad_v[1,1] - rho*grad_v[0,0]

    # Subscales definition
    vel_subscale_x = tau1*vel_residual_x
    vel_subscale_y = tau1*vel_residual_y
    pres_subscale = tau2*mas_residual

    # Compute the ASGS stabilization terms using the momentum and mass conservation residuals above
    if divide_by_rho:
        stabilization_functional = - q_gauss[0]*vel_subscale_y/y + grad_q[1]*vel_subscale_y + grad_q[0]*vel_subscale_x
    else:
        stabilization_functional = - rho*q_gauss[0]*vel_subscale_y/y + rho*grad_q[1]*vel_subscale_y + rho*grad_q[0]*vel_subscale_x
    stabilization_functional += rho*(grad_w[1,0]*v_conv_gauss[1]+w_gauss[0]*grad_v_conv[1,1])*vel_subscale_x + rho*(grad_w[0,0]*v_conv_gauss[0]+w_gauss[0]*grad_v_conv[0,0])*vel_subscale_x + grad_w[0,0]*pres_subscale
    stabilization_functional += rho*(grad_w[1,1]*v_conv_gauss[1]+w_gauss[1]*grad_v_conv[1,1])*vel_subscale_y + rho*(grad_w[0,1]*v_conv_gauss[0]+w_gauss[1]*grad_v_conv[0,0])*vel_subscale_y + grad_w[1,1]*pres_subscale + (w_gauss[1]/y)*pres_subscale
    if add_incompressibility_error:
        stabilization_functional -= 2.0*mu*w_gauss[1]*vel_subscale_y/y**2

    ## Add the stabilization terms to the original residual terms
    functional = galerkin_functional
    if ASGS_stabilization:
        functional += stabilization_functional

    ## Define DOFs and test function vectors
    dofs = sympy.zeros(n_nodes*(dim+1), 1)
    test_func = sympy.zeros(n_nodes*(dim+1), 1)

    for i in range(n_nodes):
        # Velocity DOFs and test functions
        for k in range(0,dim):
            dofs[i*(dim+1)+k] = v[i,k]
            test_func[i*(dim+1)+k] = w[i,k]
        # Pressure DOFs and test functions
        dofs[i*(dim+1)+dim] = p[i,0]
        test_func[i*(dim+1)+dim] = q[i,0]

    ## Compute LHS and RHS
    print(f"Computing {dim}D{n_nodes}N RHS Gauss point contribution\n")
    aux_functional = sympy.Matrix([[functional]])
    rhs = Compute_RHS(aux_functional.copy(), test_func, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(gauss_weight * rhs, "rRHS", mode, assignment_op='+=')

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    print(f"Computing {dim}D{n_nodes}N LHS Gauss point contribution\n")
    lhs = Compute_LHS(rhs, test_func, dofs, do_simplifications)
    lhs_out = OutputMatrix_CollectingFactors(gauss_weight * lhs, "rLHS", mode, assignment_op='+=')

    ## Replace the computed RHS and LHS in the template outstring
    outstring = outstring.replace(f"//substitute_lhs_{dim}D{n_nodes}N", lhs_out)
    outstring = outstring.replace(f"//substitute_rhs_{dim}D{n_nodes}N", rhs_out)

## Write the modified template
print("Writing output file \'" + output_filename + "\'")
out = open(output_filename,'w')
out.write(outstring)
out.close()
