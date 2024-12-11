import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Symbolic generation settings
mode = "c"
# divide_by_rho = True                # Divide the mass conservation equation by rho
ASGS_stabilization = True           # Consider ASGS stabilization terms
# add_incompressibility_error = False     # Add the incompressibility error to the viscous stress response

do_simplifications = False
output_filename = "low_mach_navier_stokes.cpp"
template_filename = "low_mach_navier_stokes_cpp_template.cpp"

info_msg = f"\n"
info_msg += f"Element generator settings:\n"
info_msg += f"\t - ASGS stabilization: {ASGS_stabilization}\n"
# info_msg += f"\t - Divide mass conservation by rho: {divide_by_rho}\n"
# info_msg += f"\t - Add incompressibility error to viscous stress: {add_incompressibility_error}\n"
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
    gauss_weight = sympy.Symbol('gauss_weight', positive = True) # Integration weight defined at cpp
    N,DN = DefineShapeFunctions(n_nodes, dim, impose_partion_of_unity, shape_functions_name='r_N', first_derivatives_name='r_DN')

    ## Unknown fields definition
    u = DefineMatrix('r_u',n_nodes,dim)            # Current step velocity (v(i,j) refers to velocity of node i component j)
    u_n = DefineMatrix('r_u_n',n_nodes,dim)        # Previous step velocity
    u_nn = DefineMatrix('r_u_nn',n_nodes,dim)      # 2 previous step velocity

    p = DefineVector('r_p',n_nodes)                # Pressure

    t = DefineVector('r_t',n_nodes)                # Temperature
    t_n = DefineVector('r_t_n',n_nodes)            # Previous step temperature
    t_nn = DefineVector('r_t_nn',n_nodes)          # 2 previous step temperature

    ## Thermodynamic pressure definition (constant in the entire domain and computed by the resolution strategy)
    p_th = sympy.Symbol('p_th', positive = True)         # Thermodynamic pressure
    dp_th_dt = sympy.Symbol('dp_th_dt', positive = True) # Thermodynamic pressure time derivative

    ## Fluid properties
    mu = sympy.Symbol('mu', positive = True)       # Dynamic viscosity
    c_p = sympy.Symbol('c_p', positive = True)     # Specific heat at constant pressure
    kappa = sympy.Symbol('kappa', positive = True) # Thermal conductivity
    alpha = sympy.Symbol('alpha', positive = True) # Thermal expansion coefficient
    rho = DefineVector('r_rho',n_nodes)              # Density
    rho_n = DefineVector('r_rho_n',n_nodes)            # Previous step density
    rho_nn = DefineVector('r_rho_nn',n_nodes)           # 2 previous step density

    ## Test functions definition
    v = DefineMatrix('v',n_nodes, dim)             # Velocity field test function
    q = DefineVector('q',n_nodes)                  # Pressure field test function
    w = DefineVector('w',n_nodes)                  # Temperature field test function

    ## Material response symbols definition
    strain_size = 3 if dim == 2 else 6
    C = DefineSymmetricMatrix('r_C', strain_size, strain_size)
    stress = DefineVector('r_stress', strain_size)

    ## Other data definitions
    g = DefineMatrix('r_g',n_nodes,dim)              # Gravity (velocity volume forcing term)
    heat_fl = DefineVector('r_heat_fl', n_nodes) # Heat flux (temperature volume forcing term)

    ## Other simbols definition
    h = sympy.Symbol('h', positive = True)      # Element size
    dt  = sympy.Symbol('dt', positive = True)   # Time increment
    dyn_tau = sympy.Symbol('dyn_tau', positive = True)
    stab_c1 = sympy.Symbol('stab_c1', positive = True)
    stab_c2 = sympy.Symbol('stab_c2', positive = True)
    stab_c3 = sympy.Symbol('stab_c3', positive = True)

    ## Backward differences coefficients
    bdf0 = sympy.Symbol('bdf0')
    bdf1 = sympy.Symbol('bdf1')
    bdf2 = sympy.Symbol('bdf2')

    ## Data interpolation to the Gauss points
    # Note that in this interpolation we do not interpolate the momentum (rho*v product) as we do not expect shocks in here
    # In other words, doing the product of the interpolations is equivalent to the interpolation of the product
    u_gauss = u.transpose()*N
    p_gauss = p.transpose()*N
    t_gauss = t.transpose()*N

    v_gauss = v.transpose()*N
    q_gauss = q.transpose()*N
    w_gauss = w.transpose()*N

    g_gauss = g.transpose()*N
    heat_fl_gauss = heat_fl.transpose()*N

    rho_gauss = (rho.transpose()*N)[0]

    ## Convective velocity definition
    #TODO: Note that we can include the linearised subscale in here to make it "more LES"
    u_conv = DefineMatrix('u_conv', n_nodes, dim) # Convective velocity defined a symbol (this is required to avoid its differentiation)
    u_conv_gauss = u_conv.transpose()*N

    ## Compute the stabilization parameters
    stab_norm_a = 0.0
    for i in range(0, dim):
        stab_norm_a += u_conv_gauss[i]**2
    stab_norm_a = sympy.sqrt(stab_norm_a)
    tau_c = 0.0 # Stabilization parameter 1
    tau_u = 0.0                                  # Stabilization parameter 2
    tau_t = 0.0                                  # Stabilization parameter 3
    # tau1 = 1.0/((rho*dyn_tau)/dt + (stab_c2*rho*stab_norm_a)/h + (stab_c1*mu)/(h*h)) # Stabilization parameter 1
    # tau2 = mu + (stab_c2*rho*stab_norm_a*h)/stab_c1                                  # Stabilization parameter 2
    # tau3 = mu + (stab_c2*rho*stab_norm_a*h)/stab_c1                                  # Stabilization parameter 3

    ## Compute the rest of magnitudes at the Gauss points
    drho_dt_gauss = ((bdf0*rho + bdf1*rho_n + bdf2*rho_nn).transpose()*N)
    du_dt_gauss = (bdf0*u + bdf1*u_n + bdf2*u_nn).transpose()*N
    dt_dt_gauss = (bdf0*t + bdf1*t_n + bdf2*t_nn).transpose()*N

    ## Gradients computation (fluid dynamics gradient)
    grad_u = DfjDxi(DN,u)
    grad_p = DfjDxi(DN,p)
    grad_t = DfjDxi(DN,t)

    grad_q = DfjDxi(DN,q)
    grad_v = DfjDxi(DN,v)
    grad_w = DfjDxi(DN,w)

    grad_rho = DfjDxi(DN,rho)

    ## Symmetric gradients calculation (already in Voigt notation)
    grad_sym_u = grad_sym_voigtform(DN,u) # Symmetric gradient of u in Voigt notation
    grad_sym_v = grad_sym_voigtform(DN,v) # Symmetric gradient of v in Voigt notation

    ## Divergences calculation
    div_u = div(DN,u)
    div_v = div(DN,v)

    ## Navier-Stokes functional
    # Mass conservation residual
    galerkin_functional = -q_gauss * drho_dt_gauss
    galerkin_functional += grad_q.transpose() * (rho_gauss * u_gauss)

    # Momentum conservation residual
    conv_term_u_gauss = u_conv_gauss.transpose() * grad_u
    galerkin_functional += v_gauss.transpose() * (rho_gauss * g_gauss)
    galerkin_functional -= v_gauss.transpose() * (rho_gauss * du_dt_gauss)
    galerkin_functional -= v_gauss.transpose() * (rho_gauss * conv_term_u_gauss.transpose())
    galerkin_functional -= grad_sym_v.transpose() * stress
    galerkin_functional += div_v * p_gauss

    # Energy conservation residual
    conv_term_t_gauss = u_conv_gauss.transpose() * grad_t
    galerkin_functional += w_gauss * heat_fl_gauss
    galerkin_functional -= w_gauss * (rho_gauss * c_p * dt_dt_gauss)
    galerkin_functional -= w_gauss * (rho_gauss * c_p * conv_term_t_gauss)
    galerkin_functional -= w_gauss * kappa * grad_w.transpose() * grad_t
    galerkin_functional += w_gauss * alpha * t_gauss * dp_th_dt

    ## Subscales definition
    R_c = - drho_dt_gauss - rho_gauss * div_u - grad_rho.transpose() * u_gauss
    R_u = rho_gauss * (g_gauss - du_dt_gauss - conv_term_u_gauss.transpose()) - grad_p
    R_t = heat_fl_gauss - rho_gauss * c_p * (dt_dt_gauss + conv_term_t_gauss) + alpha * t_gauss * dp_th_dt

    p_subscale = tau_c * R_c
    u_subscale = tau_u * R_u
    t_subscale = tau_t * R_t

    ##  Stabilization functional terms
    # Mass conservation residual
    stabilization_functional = grad_q.transpose() * (rho_gauss * u_subscale)

    # Momentum conservation residual
    conv_term_v_gauss = u_conv_gauss.transpose() * grad_v
    stabilization_functional -= drho_dt_gauss * v_gauss.transpose() * u_subscale
    stabilization_functional += rho_gauss * conv_term_v_gauss * u_subscale
    stabilization_functional += div_v * p_subscale

    # Energy conservation residual
    conv_term_w_gauss = u_conv_gauss.transpose() * grad_w
    stabilization_functional += w_gauss * alpha * t_subscale * dp_th_dt
    stabilization_functional -= c_p * drho_dt_gauss * w_gauss * t_subscale
    stabilization_functional += rho_gauss * c_p * conv_term_w_gauss * t_subscale

    ## Add the stabilization terms to the original residual terms
    functional = galerkin_functional
    if ASGS_stabilization:
        functional += stabilization_functional

    ## Multiply the functional by the corresponding integration point weight
    functional *= gauss_weight

    ## Define DOFs and test function vectors
    dofs = sympy.zeros(n_nodes*(dim+2), 1)
    test_func = sympy.zeros(n_nodes*(dim+2), 1)

    block_size = dim + 2
    for i in range(n_nodes):
        # Pressure DOFs and test functions
        p_row = i * block_size
        dofs[p_row] = p[i]
        test_func[p_row] = q[i]

        # Velocity DOFs and test functions
        u_row = i * block_size + 1
        for k in range(dim):
            dofs[u_row + k] = u[i,k]
            test_func[u_row + k] = v[i,k]
        # Temperature DOFs and test functions
        t_row = i * block_size + (dim + 1)
        dofs[t_row] = p[i,0]
        test_func[t_row] = q[i,0]

    ## Compute LHS and RHS
    print(f"Computing {dim}D{n_nodes}N RHS Gauss point contribution\n")
    aux_functional = sympy.Matrix([[functional]])
    rhs = Compute_RHS(aux_functional.copy(), test_func, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(gauss_weight * rhs, "rRHS", mode, assignment_op='+=')

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    print(f"Computing {dim}D{n_nodes}N LHS Gauss point contribution\n")
    SubstituteMatrixValue(rhs, stress, C*grad_sym_u)
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
