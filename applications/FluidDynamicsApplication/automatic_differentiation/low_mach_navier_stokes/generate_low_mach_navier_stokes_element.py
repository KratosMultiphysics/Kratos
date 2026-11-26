import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Symbolic generation settings
mode = "c"
ASGS_stabilization = True # Consider ASGS stabilization terms

do_simplifications = False
output_filename = "low_mach_navier_stokes.cpp"
template_filename = "low_mach_navier_stokes_cpp_template.cpp"

info_msg = f"\n"
info_msg += f"Element generator settings:\n"
info_msg += f"\t - ASGS stabilization: {ASGS_stabilization}\n"
print(info_msg)

dim_vector = [2, 2, 3, 3]
n_nodes_vector = [3, 4, 4, 8] # tria, quad, tet, hexa

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
    p_n = DefineVector('r_p_N',n_nodes)            # Previous step pressure
    p_nn = DefineVector('r_p_nn',n_nodes)          # 2 previous step pressure

    t = DefineVector('r_t',n_nodes)                # Temperature
    t_n = DefineVector('r_t_n',n_nodes)            # Previous step temperature
    t_nn = DefineVector('r_t_nn',n_nodes)          # 2 previous step temperature

    ## Thermodynamic pressure definition (constant in the entire domain and computed by the resolution strategy)
    p_th = sympy.Symbol('p_th', positive = True)         # Thermodynamic pressure
    dp_th_dt = sympy.Symbol('dp_th_dt', positive = True) # Thermodynamic pressure time derivative

    ## Fluid properties
    c_p = sympy.Symbol('c_p', positive = True)     # Specific heat at constant pressure
    gamma = sympy.Symbol('gamma', positive = True) # Heat capacity ratio
    kappa = sympy.Symbol('kappa', positive = True) # Thermal conductivity
    sigma = sympy.Symbol('sigma', positive = True) # Resistance (permeability inverse, 1/m^2)

    ## Test functions definition
    v = DefineMatrix('v',n_nodes, dim)             # Velocity field test function
    q = DefineVector('q',n_nodes)                  # Pressure field test function
    w = DefineVector('w',n_nodes)                  # Temperature field test function

    ## Material response symbols definition
    strain_size = 3 if dim == 2 else 6
    C = DefineSymmetricMatrix('r_C', strain_size, strain_size)
    stress = DefineVector('r_stress', strain_size)

    ## Other data definitions
    g = DefineMatrix('r_g', n_nodes, dim) # Gravity (velocity volume forcing term)
    heat_fl = DefineVector('r_heat_fl', n_nodes) # Heat flux (temperature volume forcing term)
    u_sol_frac = DefineMatrix('r_u_sol_frac', n_nodes, dim) # Solid fraction velocity

    ## Stabilization operators defined as a symbol
    tau_c = sympy.Symbol('tau_c', positive = True)
    tau_u = sympy.Symbol('tau_u', positive = True)
    tau_t = sympy.Symbol('tau_t', positive = True)

    ## Backward differences coefficients
    bdf0 = sympy.Symbol('bdf0')
    bdf1 = sympy.Symbol('bdf1')
    bdf2 = sympy.Symbol('bdf2')

    ## Symbols for linearised variables
    t_lin = DefineVector('r_t_lin', n_nodes) # Temperature defined as a symbol (to avoid differentiation in the equation of state)
    rho_lin = DefineVector('rho_lin', n_nodes) # Density defined as a symbol (to avoid differentiation in the stabilization taus)
    lin_u_conv = DefineMatrix('lin_u_conv', n_nodes, dim) # Linearised convective velocity defined as a symbol (to avoid its differentiation in the convective terms)

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
    u_sol_frac_gauss = u_sol_frac.transpose()*N

    ## Convective velocity definition
    u_m = DefineMatrix('r_u_mesh',n_nodes,dim)     # Mesh velocity
    u_conv_gauss = (u - u_m).transpose()*N      # Convective velocity
    lin_u_conv_gauss = lin_u_conv.transpose()*N # Linearised convective velocity

    ## Compute the rest of magnitudes at the Gauss points
    dp_dt_gauss = (bdf0*p + bdf1*p_n + bdf2*p_nn).transpose()*N
    du_dt_gauss = (bdf0*u + bdf1*u_n + bdf2*u_nn).transpose()*N
    dt_dt_gauss = (bdf0*t + bdf1*t_n + bdf2*t_nn).transpose()*N

    ## Gradients computation (fluid dynamics gradient)
    grad_u = DfjDxi(DN,u)
    grad_p = DfjDxi(DN,p)
    grad_t = DfjDxi(DN,t)

    grad_q = DfjDxi(DN,q)
    grad_v = DfjDxi(DN,v)
    grad_w = DfjDxi(DN,w)

    ## Symmetric gradients calculation (already in Voigt notation)
    grad_sym_u = grad_sym_voigtform(DN,u) # Symmetric gradient of u in Voigt notation
    grad_sym_v = grad_sym_voigtform(DN,v) # Symmetric gradient of v in Voigt notation

    ## Divergences calculation
    div_u = div(DN,u)
    div_v = div(DN,v)
    div_u_conv = div(DN, u - u_m)
    div_lin_u_conv = div(DN,lin_u_conv)

    ## Define state equation (ideal gases) values at Gauss point
    t_lin_gauss = t_lin.transpose()*N
    grad_t_lin = DfjDxi(DN, t_lin)

    aux_state_eq = gamma / (c_p * (gamma - 1.0))
    rho_gauss = aux_state_eq * p_th * (t_lin_gauss**-1)
    grad_rho = - grad_t_lin * (aux_state_eq * p_th * (t_lin_gauss**-2))
    drho_dt_gauss = aux_state_eq * (dp_th_dt * (t_lin_gauss**-1) - p_th * dt_dt_gauss[0] * (t_lin_gauss**-2))

    alpha = t_lin_gauss**-1

    ## Navier-Stokes functional
    # Mass conservation residual
    galerkin_functional = -q_gauss * drho_dt_gauss
    galerkin_functional -= q_gauss * (grad_rho.transpose() * u_conv_gauss)
    galerkin_functional -= q_gauss * rho_gauss * div_u

    # Momentum conservation residual
    lin_conv_term_u_gauss = lin_u_conv_gauss.transpose() * grad_u
    galerkin_functional += rho_gauss * v_gauss.transpose() * g_gauss
    galerkin_functional -= rho_gauss * v_gauss.transpose() * du_dt_gauss
    galerkin_functional -= rho_gauss * v_gauss.transpose() * lin_conv_term_u_gauss.transpose()
    galerkin_functional -= grad_sym_v.transpose() * stress
    galerkin_functional += div_v * p_gauss
    galerkin_functional -= sigma * v_gauss.transpose() * (u_gauss - u_sol_frac_gauss)

    # Energy conservation residual
    conv_term_t_gauss = u_conv_gauss.transpose() * grad_t
    galerkin_functional += w_gauss * heat_fl_gauss
    galerkin_functional -= w_gauss * (rho_gauss * c_p * dt_dt_gauss)
    galerkin_functional -= w_gauss * (rho_gauss * c_p * conv_term_t_gauss)
    galerkin_functional -= grad_w.transpose() * (kappa * grad_t)
    galerkin_functional += w_gauss * alpha * t_gauss * dp_th_dt

    ## Subscales definition
    lin_conv_term_t_gauss = lin_u_conv_gauss.transpose() * grad_t
    # R_c = - drho_dt_gauss - rho_gauss * div_u - grad_rho.transpose() * u_conv_gauss # "Standard" form (div(rho·u) = rho·div(u) + grad(rho)·u)
    R_c = - drho_dt_gauss - rho_gauss * div_u + rho_gauss * alpha * grad_t.transpose() * u_conv_gauss # "Alternative" form (div(rho·u) = rho·div(u) + grad(rho)·u = rho·div(u) + grad(rho_0 - alpha * T)·u = rho·div(u) - alpha·grad(T)·u)
    R_u = (g_gauss - du_dt_gauss - lin_conv_term_u_gauss.transpose()) * rho_gauss - grad_p - sigma * (u_gauss - u_sol_frac_gauss)
    R_t = heat_fl_gauss - rho_gauss * c_p * (dt_dt_gauss + lin_conv_term_t_gauss) + alpha * t_gauss * dp_th_dt

    p_subscale = tau_c * R_c
    u_subscale = tau_u * R_u
    t_subscale = tau_t * R_t

    ## Stabilization functional terms
    # Mass conservation residual
    stabilization_functional = rho_gauss * grad_q.transpose() * u_subscale

    # Momentum conservation residual
    lin_conv_term_v_gauss = lin_u_conv_gauss.transpose() * grad_v
    stabilization_functional += (grad_rho.transpose() * lin_u_conv_gauss) * v_gauss.transpose() * u_subscale
    stabilization_functional += (rho_gauss * div_lin_u_conv) * v_gauss.transpose() * u_subscale
    stabilization_functional += rho_gauss * lin_conv_term_v_gauss * u_subscale
    stabilization_functional += div_v * p_subscale
    stabilization_functional -= sigma * v_gauss.transpose() * u_subscale

    # Energy conservation residual
    stabilization_functional += rho_gauss * c_p * grad_w.transpose() * lin_u_conv_gauss * t_subscale
    stabilization_functional += w_gauss * c_p * (grad_rho.transpose() * lin_u_conv_gauss) * t_subscale
    stabilization_functional += w_gauss * rho_gauss * c_p * div_lin_u_conv * t_subscale
    stabilization_functional += w_gauss * alpha * t_subscale * dp_th_dt

    ## Add the stabilization terms to the original residual terms
    functional = galerkin_functional
    if ASGS_stabilization:
        functional += stabilization_functional

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
        dofs[t_row] = t[i,0]
        test_func[t_row] = w[i,0]

    ## Compute LHS and RHS
    # Compute RHS (functional differenctiation w.r.t. the shape functions nodal values)
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
