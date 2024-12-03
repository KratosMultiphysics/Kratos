import sympy
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Symbolic element generation for the Topology Optimization problems for Fluid Problems
# The following script build both elements for the primal NS problem and the ADJ NS problem
# The element LHS and RHS are than saved  for both problems in the same .cpp element file

## PHYSICS EQUATIONS ASSUMPTIONS
# -DIVIDE BY RHO: the mass conservation equation is considered divided by rho

## Symbolic generation settings
mode = "c"
do_simplifications = False

# Problem Type Definitios: "NS"-> Navier-Stokes, "ADJ_NS"->ADJOINT Navier-Stokes
problem_types = ["NS", "ADJ_NS"] # vector containing the problem types

## Prepare I/O files
# Template .cpp file
template_filename = "fluid_topology_optimization_element_template.cpp" # Template file name
print("Reading template file \'"+ template_filename + "\'\n") 
templatefile = open(template_filename) # Open template file
# Outfile
outstring       = templatefile.read() # Initialize outstring
output_filename = "fluid_topology_optimization_element.cpp" # Output file name

## Problem Settings
convective_term = True # Convection or Not
stabilization   = True # Include Stabilization

# ## ARTIFICIAL COMPRESSIBILITY NOT YET IMPLEMENTED, AND POSSIBLY WILL NEVER BE
# artificial_compressibility = False # Include Artificial Compressibility

## PROBLEM DIMENSION & N° nodes per element definition
dim_vector    = [2,3]
nnodes_vector = [(dim + 1) for dim in dim_vector]

## ITERATE OVER (dim, nnodes_x_el) couples
for dim, nnodes in zip(dim_vector, nnodes_vector):

    strain_size = dim*(dim+1)//2 # define the strain matrix size

    # Constitutive matrix definition
    C = DefineSymmetricMatrix('C',strain_size,strain_size)

    # Definition of the base functions and derivatives
    impose_partion_of_unity = False
    N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

    ## Gauss weigths definition
    gauss_weight = sympy.Symbol('gauss_weight', positive = True)

    ## Material properties
    rho    = sympy.Symbol('rho', positive = True)   # Dynamic viscosity
    mu    = sympy.Symbol('mu', positive = True)     # Dynamic viscosity
    alpha = DefineVector('alpha', nnodes)           # Resistance term
    alpha_gauss = alpha.transpose()*N               # Resistance at gauss points
    alpha = alpha_gauss[0]
    
    ## |-----------------------------------------------|
    ## | START PRIMAL NAVIER-STOKES ELEMENT DEFINITION |
    print("---\n-| NAVIER-STOKES WITH BRINKMAN PENALIZATION ELEMENT\n---")

    # Unknowns definition
    v = DefineMatrix('v',nnodes,dim)    # Current step velocity (v(i,j) refers to velocity of node i component j) 
    p = DefineVector('p',nnodes)        # Pressure

    # Forcing definitions
    f = DefineMatrix('f', nnodes, dim)  # Forcing

    # Stress definition
    stress = DefineVector('stress',strain_size)

    # Test functions definition
    w = DefineMatrix('w',nnodes,dim)    # Velocity test function
    q = DefineVector('q',nnodes)        # Pressure test function

    # Compute functions at Gauss points
    v_gauss = v.transpose()*N           # Velocity at gauss points
    p_gauss = p.transpose()*N           # Pressure at gauss points
    f_gauss = f.transpose()*N           # Forcing at gauss points
    w_gauss = w.transpose()*N           # Velocity test function at gauss points
    q_gauss = q.transpose()*N           # Pressure test function at gauss points

    # Derivatives evaluation
    # Gradient
    grad_v = DfiDxj(DN,v)   # Velocity gradient
    grad_p = DfjDxi(DN,p)   # Pressure gradient VERTICAL VECTOR
    grad_w = DfiDxj(DN,w)   # Velocity test function gradient
    grad_q = DfjDxi(DN,q)   # Pressure test function gradient VERTICAL VECTOR

    # Symmetric Gradient: 1/2*(grad(u)+grad(u)^T)
    grad_sym_v = grad_sym_voigtform(DN,v)     # Symmetric velocity gradient in Voigt notation
    grad_sym_w = grad_sym_voigtform(DN,w)     # Symmetric velocity test function gradient of w in Voigt notation

    #Divergence
    div_v = div(DN,v)   # Velocity divergence
    div_w = div(DN,w)   # Velocity test function divergence

    ## BASE GALERKIN FUNCTIONAL RESIDUAL
    rv_galerkin  = -(grad_sym_w.transpose()*stress) # diffusion term, OBS--> A::B=A_sym::B_sym + A_asym::B_asym, but (grad(u)+grad(u)^T) is sym, so we can use the grad_sym(w)
    rv_galerkin +=  div_w*p_gauss # pressure term
    rv_galerkin += -alpha*(w_gauss.transpose()*v_gauss) # resistance (alpha) term
    rv_galerkin +=  rho*w_gauss.transpose()*f_gauss # RHS forcing
    rv_galerkin += -q_gauss*div_v # mass conservation term (divided by rho)
    rv = rv_galerkin # save BASE residual into TOTAL element residual

    ## HANDLE CONVECTION 
    if (convective_term):
        # Definition
        vconv = DefineMatrix('vconv',nnodes,dim)    # Convective velocity 
        # Gauss points
        vconv_gauss = vconv.transpose()*N           # Convective velocity on gauus points
        # Derivatives
        div_vconv = div(DN, vconv)                 # Divergence of the convective velocity
        # CONVECTIVE GALERKIN FUNCTIONAL RESIDUAL
        rv_conv = -rho*w_gauss.transpose()*(grad_v*vconv_gauss) # convection term
        rv += rv_conv # save CONVECTION residual into TOTAL element residual
    # END HANDLE CONVECTION

    ## HANDLE STABILIZATION
    if (stabilization):
        # Parameters definition
        h = sympy.Symbol('h', positive = True)      # Mesh size
        # dt = sympy.Symbol('dt', positive = True) # Time increment
        # dyn_tau = sympy.Symbol('dyn_tau', positive = True) # Time stabilization constant
        #Stabilization constants definition
        stab_c1 = sympy.Symbol('stab_c1', positive = True) # Diffusion stabilization constant
        stab_c2 = sympy.Symbol('stab_c2', positive = True) # Convection stabilization constant
        stab_c3 = sympy.Symbol('stab_c3', positive = True) # Resistance (alpha) term stabilization constant
        # TAU1 EVALUATION steps
        tau1_denominator = 0.0
        # tau1_denominator += rho*dyn_tau/dt # time contribution to tau_1
        tau1_denominator += stab_c1*mu/h**2 # diffusion contribution to tau_1
        tau1_denominator += stab_c3*alpha # resistance contribution to tau_1
        # TAU2 EVALUATION steps
        tau2_factor  = 0.0
        tau2_factor += stab_c3*alpha*h*h # resistance contribution to tau_2
        # Update TAU1 & TAU2 with CONVECTIVE term stabilization
        if (convective_term):
            stab_norm_a = 0.0
            for i in range(0, dim): # evaluate the vconv norm
                stab_norm_a += vconv_gauss[i]**2
            stab_norm_a = sympy.sqrt(stab_norm_a)
            tau1_denominator += stab_c2*rho*stab_norm_a/h # convection contribution to tau_1
            tau2_factor      += stab_c2*rho*stab_norm_a*h # convection contribution to tau_2
        # Definition of TAU1 & TAU2
        tau1 = 1.0/tau1_denominator         # Stabilization parameter 1
        tau2 = mu + tau2_factor/stab_c1     # Stabilization parameter 2
        ##  STABILIZATION functional terms
        # momentum conservation residual
        momentum_residual  =  rho*f_gauss # forcing
        momentum_residual += -grad_p # pressure gradient
        # diffusion eliminates since it involves 2nd order derivatives for 1st order solutions
        momentum_residual += -(alpha)*v_gauss
        if (convective_term):
            momentum_residual += -rho*grad_v*vconv_gauss
        # mass conservation residual
        mass_residual = -div_v
        # Velocity and Presure subscales
        velocity_subscale = tau1*momentum_residual
        pressure_subscale = tau2*mass_residual
        # ASGS stabilization
        # some of the terms in rv_stab are taken with sign (+) since they derive from an integration by parts 
        rv_stab  =  velocity_subscale.transpose()*grad_q
        rv_stab += -alpha*velocity_subscale.transpose()*w_gauss # stab resistance residual
        rv_stab +=  pressure_subscale*div_w
        if (convective_term):
            rv_stab += rho*velocity_subscale.transpose()*(grad_w*vconv_gauss) # 1st stab convective residual
            rv_stab += rho*div_vconv*velocity_subscale.transpose()*w_gauss  # 2nd stab convective residual
        rv += rv_stab # save STABILIZATION residual into TOTAL element residual
    # END HANDLE STABILIZATION

    ## Define DOFs and test function vectors
    dofs = sympy.zeros(nnodes*(dim+1), 1)
    testfunc = sympy.zeros(nnodes*(dim+1), 1)
    for i in range(nnodes):
        # Velocity DOFs and test functions
        for k in range(0,dim):
            dofs[i*(dim+1)+k] = v[i,k]
            testfunc[i*(dim+1)+k] = w[i,k]
        # Pressure DOFs and test functions
        dofs[i*(dim+1)+dim] = p[i,0]
        testfunc[i*(dim+1)+dim] = q[i,0] 

    # COMPUTE RHS AND LHS
    print(f"Computing {dim}D{nnodes}N NS RHS Gauss point contribution\n")
    rhs = Compute_RHS(rv.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(gauss_weight*rhs, "rRHS", mode, assignment_op='+=')
    print(f"Computing {dim}D{nnodes}N NS LHS Gauss point contribution\n")
    SubstituteMatrixValue(rhs, stress, C*grad_sym_v)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = OutputMatrix_CollectingFactors(gauss_weight*lhs, "rLHS", mode, assignment_op='+=')
    # print(OutputMatrix_CollectingFactors(lhs,"lhs",mode))
    ## Replace the computed RHS and LHS in the template outstring
    outstring = outstring.replace(f"//substitute_lhs_{dim}D{nnodes}N_NS", lhs_out)
    outstring = outstring.replace(f"//substitute_rhs_{dim}D{nnodes}N_NS", rhs_out)
    ## | END PRIMAL NAVIER-STOKES ELEMENT DEFINITION |
    ## |---------------------------------------------|
    
    ## |------------------------------------------------|
    ## | START ADJOINT NAVIER-STOKES ELEMENT DEFINITION |
    print("---\n-| ADJOINT NAVIER-STOKES WITH BRINKMAN PENALIZATION ELEMENT\n---")

    # Unknowns definition
    v_adj = DefineMatrix('v_adj',nnodes,dim)    # Current step velocity_adj (v_adj(i,j) refers to velocity_adj of node i component j) 
    p_adj = DefineVector('p_adj',nnodes)        # Pressure_adj

    # Forcing_adj definitions
    f_adj = DefineMatrix('f_adj', nnodes, dim)  # Forcing_adj

    # Stress_adj definition
    stress_adj = DefineVector('stress_adj',strain_size)

    # Test functions definition
    w_adj = DefineMatrix('w_adj',nnodes,dim)    # Velocity_adj test function
    q_adj = DefineVector('q_adj',nnodes)        # Pressure_adj test function

    # Compute functions at Gauss points
    v_adj_gauss = v_adj.transpose()*N           # Velocity_adj at gauss points
    p_adj_gauss = p_adj.transpose()*N           # Pressure_adj at gauss points
    f_adj_gauss = f_adj.transpose()*N           # Forcing_adj at gauss points
    w_adj_gauss = w_adj.transpose()*N           # Velocity_adj test function at gauss points
    q_adj_gauss = q_adj.transpose()*N           # Pressure_adj test function at gauss points
    # RESITANCE ALPHA IS THE SAME OF THE PRIMAL
    # alpha_adj_gauss = alpha.transpose()*N   # Resistance at gauss points

    # Derivatives evaluation
    # Gradient
    grad_v_adj = DfiDxj(DN,v_adj)   # Velocity_adj gradient
    grad_p_adj = DfjDxi(DN,p_adj)   # Pressure_adj gradient VERTICAL VECTOR
    grad_w_adj = DfiDxj(DN,w_adj)   # Velocity_adj test function gradient
    grad_q_adj = DfjDxi(DN,q_adj)   # Pressure_adj test function gradient VERTICAL VECTOR

    # Symmetric Gradient: 1/2*(grad(u)+grad(u)^T)
    grad_sym_v_adj = grad_sym_voigtform(DN,v_adj)     # Symmetric velocity_adj gradient in Voigt notation
    grad_sym_w_adj = grad_sym_voigtform(DN,w_adj)     # Symmetric velocity_adj test function gradient of w in Voigt notation

    #Divergence
    div_v_adj = div(DN,v_adj)   # Velocity_adj divergence
    div_w_adj = div(DN,w_adj)   # Velocity_adj test function divergence

    ## ADJOINT BASE GALERKIN FUNCTIONAL RESIDUAL
    rv_galerkin_adj  = -(grad_sym_w_adj.transpose()*stress_adj) # diffusion term, OBS--> A::B=A_sym::B_sym + A_asym::B_asym, but (grad(u)+grad(u)^T) is sym, so we can use the grad_sym(w)
    rv_galerkin_adj +=  div_w_adj*p_adj_gauss # pressure term
    rv_galerkin_adj += -alpha*(w_adj_gauss.transpose()*v_adj_gauss) # resistance (alpha) term
    rv_galerkin_adj +=  rho*w_adj_gauss.transpose()*f_adj_gauss # RHS forcing
    rv_galerkin_adj += -q_adj_gauss*div_v_adj # mass conservation term
    rv_adj = rv_galerkin_adj # save BASE residual into TOTAL element residual

    ## HANDLE ADJOINT CONVECTION 
    if (convective_term):
        # Definition
        vconv_adj = DefineMatrix('vconv_adj',nnodes,dim)    # Convective_adj velocity 
        # Gauss points
        vconv_adj_gauss = vconv_adj.transpose()*N       # Convective_adj velocity on gauus points
        # Derivatives
        grad_vconv_adj = DfiDxj(DN,vconv_adj)           # Gradient of the convective_adj velocity
        div_vconv_adj  = div(DN, vconv_adj)              # Divergence of the convective_adj velocity
        # CONVECTIVE ADJOINT GALERKIN FUNCTIONAL RESIDUAL
        # Adjoint velocity convection driven by the vconv = NS velocity 
        # -> evaluated in a weak form such that the Neumann BC are imposed "naturally"
        # from the Adjoint PDE its formulation should be: rv_galerkin  += rho*(w_gauss.transpose()*(grad_v*vconv_gauss))
        rv_conv_adj   = -rho*(v_adj_gauss.transpose()*(grad_w_adj*vconv_adj_gauss)) 
        # Velocity convection driven by the adjoint velocity: called ADJOINT RESISTANCE caouse we can write (alphaI+grad(U))u_a * w
        # additional term appearing in the ADJ_NS system w.r.t. the NS eqs.
        rv_conv_adj  += -rho*(w_adj_gauss.transpose()*(grad_vconv_adj*v_adj_gauss))
        rv_adj += rv_conv_adj # save CONVECTION residual into TOTAL element residual
    # END HANDLE CONVECTION

    # HANDLE ADJOINT STABILIZATION
    if (stabilization):
        # Parameters definition
        h = sympy.Symbol('h', positive = True)      # Mesh size
        # dt = sympy.Symbol('dt', positive = True) # Time increment
        # dyn_tau = sympy.Symbol('dyn_tau', positive = True) # Time stabilization constant
        #Stabilization constants definition
        stab_c1 = sympy.Symbol('stab_c1', positive = True) # Diffusion stabilization constant
        stab_c2 = sympy.Symbol('stab_c2', positive = True) # Convection stabilization constant
        stab_c3 = sympy.Symbol('stab_c3', positive = True) # Resistance (alphaI+grad(u)) term stabilization constant
        # TAU1 EVALUATION steps
        tau1_denominator_adj = 0.0
        # tau1_denominator += rho*dyn_tau/dt # time contribution to tau_1
        tau1_denominator_adj += stab_c1*mu/h**2 # diffusion contribution to tau_1
        tau1_denominator_adj += stab_c3*alpha # resistance (alpha) contribution to tau_1
        # TAU2 EVALUATION steps
        tau2_factor_adj  = 0.0
        tau2_factor_adj += stab_c3*alpha*h*h # resistance contribution to tau_2
        # Update TAU1 & TAU2 with CONVECTIVE term stabilization
        if (convective_term):
            # USED:     the true convection could be stabilized as the classic convection term with ||a||
            #           the term similar to the resistance term stabilization has been added to the Resistance stabilization using the frobenius norm of ||grad(U)||_frob
            # TRUE CONVECTION STABILIZATION
            stab_norm_a_adj = 0.0
            for i in range(0, dim): # evaluate the vconv norm
                stab_norm_a_adj += vconv_adj_gauss[i]**2
            stab_norm_a_adj = sympy.sqrt(stab_norm_a_adj)
            tau1_denominator_adj += stab_c2*rho*stab_norm_a_adj/h # convection contribution to tau_1
            tau2_factor_adj += stab_c2*rho*stab_norm_a_adj*h # convection contribution to tau_2
            # grad(U) ADJOINT RESISTANCE STABILIZATION
            stab_norm_gradU_adj = 0.0
            for i in range(0,dim):
                for j in range(0,dim):
                    stab_norm_gradU_adj += grad_vconv_adj[i,j]**2
            stab_norm_gradU_adj = sympy.sqrt(stab_norm_a_adj)
            tau1_denominator_adj += stab_c3*stab_norm_gradU_adj # convection contribution to tau_1
            tau2_factor_adj += stab_c3*stab_norm_gradU_adj*h*h # convection contribution to tau_2
        # Definition of TAU1 & TAU2
        tau1_adj = 1.0/tau1_denominator_adj         # Stabilization parameter 1
        tau2_adj = mu + tau2_factor_adj/stab_c1     # Stabilization parameter 2
        ##  STABILIZATION functional terms
        # momentum conservation residual
        momentum_residual_adj  =  rho*f_adj_gauss # forcing
        momentum_residual_adj += -grad_p_adj # pressure gradient
        # diffusion eliminates since it involves 2nd order derivatives for 1st order solutions
        momentum_residual_adj += -alpha*v_adj_gauss
        if (convective_term):
            momentum_residual_adj +=  rho*grad_v_adj*vconv_adj_gauss
            momentum_residual_adj += -rho*grad_vconv_adj*v_adj_gauss
        # mass conservation residual
        mass_residual_adj = -div_v_adj
        # Velocity and Presure subscales
        velocity_adj_subscale = tau1_adj*momentum_residual_adj
        pressure_adj_subscale = tau2_adj*mass_residual_adj
        # ASGS stabilization
        # some of the terms in rv_stab are taken with sign (+) since they derive from an integration by parts 
        rv_stab_adj  =  velocity_adj_subscale.transpose()*grad_q_adj
        rv_stab_adj += -alpha*velocity_adj_subscale.transpose()*w_adj_gauss # stab resistance residual
        rv_stab_adj +=  pressure_adj_subscale*div_w_adj
        if (convective_term):
            rv_stab_adj += -rho*velocity_adj_subscale.transpose()*(grad_w_adj*vconv_adj_gauss) # 1st stab convective residual: convective term
            rv_stab_adj += -rho*w_adj_gauss.transpose()*(grad_vconv_adj*velocity_adj_subscale)  # 2nd stab convective residual: adjoint resistance term
        rv_adj += rv_stab_adj # save STABILIZATION residual into TOTAL element residual
    # END HANDLE STABILIZATION

    ## Define DOFs and test function vectors
    dofs = sympy.zeros(nnodes*(dim+1), 1)
    testfunc = sympy.zeros(nnodes*(dim+1), 1)
    for i in range(nnodes):
        # Velocity DOFs and test functions
        for k in range(0,dim):
            dofs[i*(dim+1)+k] = v_adj[i,k]
            testfunc[i*(dim+1)+k] = w_adj[i,k]
        # Pressure DOFs and test functions
        dofs[i*(dim+1)+dim] = p_adj[i,0]
        testfunc[i*(dim+1)+dim] = q_adj[i,0] 

    # COMPUTE RHS AND LHS
    print(f"Computing {dim}D{nnodes}N ADJ-NS RHS Gauss point contribution\n")
    rhs_adj = Compute_RHS(rv_adj.copy(), testfunc, do_simplifications)
    rhs_adj_out = OutputVector_CollectingFactors(gauss_weight*rhs_adj, "rRHS", mode, assignment_op='+=')
    print(f"Computing {dim}D{nnodes}N ADJ-NS LHS Gauss point contribution\n")
    SubstituteMatrixValue(rhs_adj, stress_adj, C*grad_sym_v_adj)
    lhs_adj = Compute_LHS(rhs_adj, testfunc, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_adj_out = OutputMatrix_CollectingFactors(gauss_weight*lhs_adj, "rLHS", mode, assignment_op='+=')
    # print(OutputMatrix_CollectingFactors(lhs,"lhs",mode))
    ## Replace the computed RHS and LHS in the template outstring
    outstring = outstring.replace(f"//substitute_lhs_{dim}D{nnodes}N_ADJ_NS", lhs_adj_out)
    outstring = outstring.replace(f"//substitute_rhs_{dim}D{nnodes}N_ADJ_NS", rhs_adj_out)
    ## | END ADJOINT NAVIER-STOKES ELEMENT DEFINITION |
    ## |----------------------------------------------|

## Write the modified template
print("Writing output file \'" + output_filename + "\'")
out = open(output_filename,'w')
out.write(outstring)
out.close()

