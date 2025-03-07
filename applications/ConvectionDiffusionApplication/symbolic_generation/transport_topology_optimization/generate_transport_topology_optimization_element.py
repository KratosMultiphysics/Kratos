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
problem_types = ["T", "ADJ_T"] # vector containing the problem types

## Prepare I/O files
# Template .cpp file
template_filename = "transport_topology_optimization_element_template.cpp" # Template file name
print("Reading template file \'"+ template_filename + "\'\n") 
templatefile = open(template_filename) # Open template file
# Outfile
outstring       = templatefile.read() # Initialize outstring
output_filename = "transport_topology_optimization_element.cpp" # Output file name

## Problem Settings
convective_term = True # Convective Term
stabilization   = True # Include Stabilization
include_functionals = True # Include Functional Terms in the Adjoint Source Term Term
n_functionals = 20 # Number of implemented functionals

## PROBLEM DIMENSION & N° nodes per element definition
# dim_vector    = [2,3]
# nnodes_vector = [(dim + 1) for dim in dim_vector]
dim_vector = [2, 2, 3]
nnodes_vector = [3, 4, 4] # tria, quad, tet

## ITERATE OVER (dim, nnodes_x_el) couples
for dim, nnodes in zip(dim_vector, nnodes_vector):

    # Definition of the base functions and derivatives
    impose_partion_of_unity = False
    N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

    ## Gauss weigths definition
    gauss_weight = sympy.Symbol('gauss_weight', positive = True)

    ## Material properties
    Conductivity       = DefineVector('D',nnodes)    # Conductivity
    Conductivity_gauss = Conductivity.transpose()*N  # Conductivity at gauss points

    Decay        = DefineVector('k',nnodes)    # Decay Factor
    Decay_gauss  = Decay.transpose()*N  # Decay Factor at gauss points

    ConvCoeff        = DefineVector('c',nnodes)    # convection coefficient
    ConvCoeff_gauss  = ConvCoeff.transpose()*N  # convection coefficient at gauss points

    # print(ConvCoeff_gauss.shape)
    # print(ConvCoeff.shape)
    # input()
    # decay_gauss * (t_gauss.transpose()*q_gauss)
    
    ## |-----------------------------------------------|
    ## | START PRIMAL TOP. OPT. TRANSPORT ELEMENT DEFINITION |
    print("---\n-| FLUID TOPOLOGY OPTIMIZATION TRANSPORT ELEMENT\n---")

    # Unknowns definition
    t = DefineVector('t',nnodes)        # Temperature

    # Source Term definitions
    f = DefineVector('source', nnodes)  # Source Term

    # Test functions definition
    q = DefineVector('q',nnodes)        # Temperature test function

    # Compute functions at Gauss points
    t_gauss = t.transpose()*N           # Temperature at gauss points
    f_gauss = f.transpose()*N           # Source Term at gauss points
    q_gauss = q.transpose()*N           # Temperature test function at gauss points

    # Derivatives evaluation
    # Gradient
    grad_t = DfjDxi(DN,t)   # Temperature gradient VERTICAL VECTOR
    grad_q = DfjDxi(DN,q)   # Temperature test function gradient VERTICAL VECTOR

    ## BASE GALERKIN FUNCTIONAL RESIDUAL
    rv_conductivity  = -Conductivity_gauss*(grad_q.transpose()*grad_t)  # diffusion term
    rv_source =  q_gauss.transpose()*f_gauss # RHS Source Term
    rv_decay = -Decay_gauss*(q_gauss.transpose()*t_gauss) # reaction term
    rv = rv_conductivity + rv_source + rv_decay # save BASE residual into TOTAL element residual

    ## HANDLE CONVECTION 
    if (convective_term):
        # Definition
        vconv = DefineMatrix('vconv',nnodes,dim)    # Convective velocity 
        # Gauss points
        vconv_gauss = vconv.transpose()*N           # Convective velocity on gauss points
        # CONVECTIVE GALERKIN FUNCTIONAL RESIDUAL
        rv_conv = -ConvCoeff_gauss*(q_gauss.transpose()*(vconv_gauss.transpose()*grad_t)) # convective term
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
        tau1_denominator += stab_c1*Conductivity_gauss[0]/h**2 # diffusion contribution to tau_1
        tau1_denominator += stab_c3*Decay_gauss[0] # decay contribution to tau_1
        # Update TAU1 with CONVECTIVE term stabilization
        if (convective_term):
            stab_norm_a = 0.0
            for i in range(0, dim): # evaluate the vconv norm
                stab_norm_a += vconv_gauss[i]**2
            stab_norm_a = sympy.sqrt(stab_norm_a)
            tau1_denominator += stab_c2*abs(ConvCoeff_gauss[0])*stab_norm_a/h # convection contribution to tau_1
        # Definition of TAU1 & TAU2
        tau1 = 1.0/tau1_denominator         # Stabilization parameter 1
        ##  STABILIZATION functional terms
        # mass conservation residual
        mass_residual  = f_gauss # Source Term
        mass_residual += -Decay_gauss*t_gauss # Decay_gauss factor
        # eliminates diffusion involving 2nd order derivatives for 1st order solutions
        if (convective_term):
            mass_residual += -ConvCoeff_gauss*vconv_gauss.transpose()*grad_t # convective term
        # Temperature subscales
        temperature_subscale = tau1*mass_residual
        # ASGS stabilization
        # some of the terms in rv_stab are taken with sign (+) since they derive from an integration by parts 
        rv_stab = -Decay_gauss*q_gauss.transpose()*temperature_subscale # stab decay residual
        if (convective_term):
            rv_stab += ConvCoeff_gauss*temperature_subscale*vconv_gauss.transpose()*grad_q # stab convection
        rv += rv_stab # save STABILIZATION residual into TOTAL element residual
    # END HANDLE STABILIZATION

    ## Define DOFs and test function vectors
    dofs = sympy.zeros(nnodes, 1)
    testfunc = sympy.zeros(nnodes, 1)
    for i in range(nnodes):
        # Temperature DOFs and test functions
        dofs[i] = t[i,0]
        testfunc[i] = q[i,0]

    # COMPUTE RHS AND LHS
    print(f"Computing {dim}D{nnodes}N T RHS Gauss point contribution\n")
    rhs = Compute_RHS(rv.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(gauss_weight*rhs, "rRHS", mode, assignment_op='+=')
    print(f"Computing {dim}D{nnodes}N T LHS Gauss point contribution\n")
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = OutputMatrix_CollectingFactors(gauss_weight*lhs, "rLHS", mode, assignment_op='+=')
    # print(OutputMatrix_CollectingFactors(lhs,"lhs",mode))
    ## Replace the computed RHS and LHS in the template outstring
    outstring = outstring.replace(f"//substitute_lhs_{dim}D{nnodes}N_T", lhs_out)
    outstring = outstring.replace(f"//substitute_rhs_{dim}D{nnodes}N_T", rhs_out)
    ## | END PRIMAL TOP. OPT. TRANSPORT ELEMENT DEFINITION |
    ## |---------------------------------------------|
    
    ## |------------------------------------------------|
    ## | START ADJOINT TOP. OPT. TRANSPORT ELEMENT DEFINITION |
    print("---\n-| ADJOINT FLUID TOPOLOGY OPTIMIZATION TRANSPORT ELEMENT\n---")

    # Unknowns definition
    t_adj = DefineVector('t_adj',nnodes)        # Temperature

    # Source Term definitions
    f_adj = DefineVector('source_adj', nnodes)  # Source Term

    # Test functions definition
    q_adj = DefineVector('q_adj', nnodes)        # Temperature test function

    # Compute functions at Gauss points
    t_adj_gauss = t_adj.transpose()*N           # Temperature at gauss points
    f_adj_gauss = f_adj.transpose()*N           # Source Term at gauss points
    q_adj_gauss = q_adj.transpose()*N           # Temperature test function at gauss points

    # Derivatives evaluation
    # Gradient
    grad_t_adj = DfjDxi(DN,t_adj)   # Temperature gradient VERTICAL VECTOR
    grad_q_adj = DfjDxi(DN,q_adj)   # Temperature test function gradient VERTICAL VECTOR

    ## BASE GALERKIN FUNCTIONAL RESIDUAL
    rv_conductivtiy_adj  = -Conductivity_gauss*(grad_q_adj.transpose()*grad_t_adj)  # diffusion term
    rv_source_adj =  q_adj_gauss.transpose()*f_adj_gauss # RHS Source Term
    rv_decay_adj = -Decay_gauss*(q_adj_gauss.transpose()*t_adj_gauss) # reaction term
    rv_adj = rv_conductivtiy_adj + rv_source_adj + rv_decay_adj # save BASE residual into TOTAL element residual
    
    ## HANDLE ADJOINT FORCING TERMS  
    if (include_functionals):
        # Functionals Database
        # 0: resistance  : int_{\Omega}{alpha*||u||^2}
        # 1: strain-rate : int_{\Omega}{2*mu*||S||^2} , with S = 1/2*(grad(u)+grad(u)^T) strain-rate tensor
        # 2: vorticity   : int_{\Omega}{2*mu*||R||^2} = int_{\Omega}{mu*||curl(u)||^2} , curl(u) = vorticity vector, R = 1/2*(grad(u)-grad(u)^T) rotation-rate tensor
        # 3: outlet_transport_scalar : int_{\Gamma_{out}}{c}
        # 4: region_transport_scalar: int_{\Omega}{c^2}
        # 5: transport_scalar_diffusion: int_{\Omega}{D\\||grad(u)||^2}
	    # 6: transport_scalar_convection: int_{\Omega}{beta*T*dot(u,grad(T))}
	    # 7: transport_scalar_decay: int_{\Omega}{kT^2}
	    # 8: transport_scalar_source: int_{\Omega}{-Q*T}
        # Source Term definitions
        functional_weights = DefineVector('functional_weights', n_functionals) # Weights of the functionals terms

        # VARIABLES
        opt_t = DefineVector('opt_t', nnodes)  
        t_physics = DefineVector('t', nnodes) 
        vconv_physics = DefineMatrix('vconv',nnodes,dim)
        f_physics = DefineVector('source', nnodes) 
        # VARIABLES AT GAUSS POINTS
        opt_t_gauss = opt_t.transpose()*N # Optimization Temperature at gauss points
        t_physics_gauss = t_physics.transpose()*N # Transport sclalar at gauss points
        vconv_physics_gauss = vconv_physics.transpose()*N # Physics convective velocity at gauss points
        f_physics_gauss = f_physics.transpose()*N # Physics source at gauss points
        # GRADIENTS
        grad_t_physics = DfjDxi(DN,t_physics)
        
        # region_transport_scalar functional
        rv_funct_region_transport_scalar  = 2.0*(q_adj_gauss.transpose()*opt_t_gauss)
        # transport_scalar_transfer
        rv_funct_transport_scalar_diffusion  = 2.0*Conductivity_gauss*(grad_q_adj.transpose()*grad_t_physics)
        rv_funct_transport_scalar_convection = ConvCoeff_gauss*((vconv_physics_gauss.transpose()*grad_t_physics)*q_adj_gauss + t_physics_gauss*(vconv_physics_gauss.transpose()*grad_q_adj)) 
        rv_funct_transport_scalar_decay      = 2.0*Decay_gauss*(q_adj_gauss*t_physics_gauss)
        rv_funct_transport_scalar_source     = -f_physics_gauss*q_adj_gauss
        
        # SUM FUNCTIONAL CONTRIBUTIONS TO ADJ RESIDUAL
        rv_adj -= functional_weights[4]*rv_funct_region_transport_scalar 
        rv_adj -= functional_weights[5]*rv_funct_transport_scalar_diffusion
        rv_adj -= functional_weights[6]*rv_funct_transport_scalar_convection
        rv_adj -= functional_weights[7]*rv_funct_transport_scalar_decay
        rv_adj -= functional_weights[8]*rv_funct_transport_scalar_source
    
    ## HANDLE CONVECTION 
    if (convective_term):
        # Definition
        vconv_adj = DefineMatrix('vconv_adj',nnodes,dim)    # Convective velocity 
        # Gauss points
        vconv_adj_gauss = vconv_adj.transpose()*N           # Convective velocity on gauus points
        # CONVECTIVE GALERKIN FUNCTIONAL RESIDUAL
        rv_conv_adj = ConvCoeff_gauss*(q_adj_gauss.transpose()*(vconv_adj_gauss.transpose()*grad_t_adj)) # convective term
        rv_adj += rv_conv_adj # save CONVECTION residual into TOTAL element residual
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
        stab_c3 = sympy.Symbol('stab_c3', positive = True) # Decay term stabilization constant
        # TAU1 EVALUATION steps
        tau1_denominator_adj = 0.0
        # tau1_denominator += rho*dyn_tau/dt # time contribution to tau_1
        tau1_denominator_adj += stab_c1*Conductivity_gauss[0]/h**2 # diffusion contribution to tau_1
        tau1_denominator_adj += stab_c3*Decay_gauss[0] # decay contribution to tau_1
        # Update TAU1 with CONVECTIVE term stabilization
        if (convective_term):
            stab_norm_a_adj = 0.0
            for i in range(0, dim): # evaluate the vconv norm
                stab_norm_a_adj += vconv_adj_gauss[i]**2
            stab_norm_a_adj = sympy.sqrt(stab_norm_a_adj)
            tau1_denominator_adj += stab_c2*abs(ConvCoeff_gauss[0])*stab_norm_a_adj/h # convection contribution to tau_1
        # Definition of TAU1 & TAU2
        tau1_adj = 1.0/tau1_denominator_adj         # Stabilization parameter 1
        ##  STABILIZATION functional terms
        # mass conservation residual
        mass_residual_adj  = f_adj_gauss # Source Term
        mass_residual_adj += -Decay_gauss*t_adj_gauss # decay factor
        # eliminates diffusion involving 2nd order derivatives for 1st order solutions
        if (convective_term):
            mass_residual_adj += ConvCoeff_gauss*vconv_adj_gauss.transpose()*grad_t_adj # convective term

        ## Include adjoint forcing terms coming from the functional definition
        if (include_functionals):
            mass_residual_adj += -functional_weights[4]*2.0*opt_t_gauss
            # residual coming from transfer functional DIFFUSION  term involves second order derivatives. 
            # residual coming from transfer functional CONVECTION term erases when evaluated not in the weak form, only a term in grad(beta) remains
            grad_ConvCoeff = DfjDxi(DN,ConvCoeff) 
            mass_residual_adj += functional_weights[6]*((vconv_physics_gauss.transpose()*grad_ConvCoeff)*t_physics_gauss)
            mass_residual_adj += -functional_weights[7]*(2.0*Decay_gauss*t_physics_gauss)
            mass_residual_adj += functional_weights[8]*(f_physics_gauss) 

        # Temperature subscales
        temperature_subscale_adj = tau1_adj*mass_residual_adj

        # ASGS stabilization
        # some of the terms in rv_stab are taken with sign (+) since they derive from an integration by parts 
        rv_stab_adj = -Decay_gauss*q_adj_gauss.transpose()*temperature_subscale_adj # stab decay residual
        if (convective_term):
            rv_stab_adj += -ConvCoeff_gauss*temperature_subscale_adj*vconv_adj_gauss.transpose()*grad_q_adj # stab convection
        rv_adj += rv_stab_adj # save STABILIZATION residual into TOTAL element residual
    # END HANDLE STABILIZATION

    ## Define DOFs and test function vectors
    dofs = sympy.zeros(nnodes, 1)
    testfunc = sympy.zeros(nnodes, 1)
    for i in range(nnodes):
        # Temperature adj DOFs and test functions
        dofs[i] = t_adj[i,0]
        testfunc[i] = q_adj[i,0] 

    # COMPUTE RHS AND LHS
    print(f"Computing {dim}D{nnodes}N ADJ-T RHS Gauss point contribution\n")
    rhs_adj = Compute_RHS(rv_adj.copy(), testfunc, do_simplifications)
    rhs_adj_out = OutputVector_CollectingFactors(gauss_weight*rhs_adj, "rRHS", mode, assignment_op='+=')
    print(f"Computing {dim}D{nnodes}N ADJ-T LHS Gauss point contribution\n")
    lhs_adj = Compute_LHS(rhs_adj, testfunc, dofs, do_simplifications) # Compute the LHS
    lhs_adj_out = OutputMatrix_CollectingFactors(gauss_weight*lhs_adj, "rLHS", mode, assignment_op='+=')
    # print(OutputMatrix_CollectingFactors(lhs,"lhs",mode))
    ## Replace the computed RHS and LHS in the template outstring
    outstring = outstring.replace(f"//substitute_lhs_{dim}D{nnodes}N_ADJ_T", lhs_adj_out)
    outstring = outstring.replace(f"//substitute_rhs_{dim}D{nnodes}N_ADJ_T", rhs_adj_out)
    ## | END ADJOINT TOP. OPT. TRANSPORT ELEMENT DEFINITION |
    ## |----------------------------------------------|

## Write the modified template
print("Writing output file \'" + output_filename + "\'")
out = open(output_filename,'w')
out.write(outstring)
out.close()

