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

## Prepare I/O files
# Template .cpp file
template_filename = "topology_optimization_pde_filter_element_template.cpp" # Template file name
print("Reading template file \'"+ template_filename + "\'\n") 
templatefile = open(template_filename) # Open template file
# Outfile
outstring       = templatefile.read() # Initialize outstring
output_filename = "topology_optimization_pde_filter_element.cpp" # Output file name

## Problem Settings
stabilization   = True # Include Stabilization

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

    Decay        = DefineVector('r',nnodes)    # Decay Factor
    Decay_gauss  = Decay.transpose()*N  # Decay Factor at gauss points
    
    ## |-----------------------------------------------|
    ## | START ELEMENT DEFINITION |
    print("---\n-|TOPOLOGY OPTIMIZATION PDE FILTER ELEMENT\n---")

    # Unknowns definition
    f = DefineVector('filter',nnodes)        # Filtered value

    # Source Term definitions
    s = DefineVector('source', nnodes)  # Source Term = non filtered value

    # Test functions definition
    q = DefineVector('q',nnodes)        # Filter test function

    # Compute functions at Gauss points
    f_gauss = f.transpose()*N           # filtered value at gauss points
    s_gauss = s.transpose()*N           # source term at gauss points
    q_gauss = q.transpose()*N           # Filter test function at gauss points

    # Derivatives evaluation
    # Gradient
    grad_f = DfjDxi(DN,f)   # Filter gradient VERTICAL VECTOR
    grad_q = DfjDxi(DN,q)   # Filter test function gradient VERTICAL VECTOR

    ## BASE GALERKIN FUNCTIONAL RESIDUAL
    rv_conductivity  = -Conductivity_gauss*(grad_q.transpose()*grad_f)  # diffusion term
    rv_source =  q_gauss.transpose()*s_gauss # RHS Source Term
    rv_decay = -Decay_gauss*(q_gauss.transpose()*f_gauss) # reaction term
    rv = rv_conductivity + rv_source + rv_decay # save BASE residual into TOTAL element residual

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
        # Definition of TAU1 & TAU2
        tau1 = 1.0/tau1_denominator         # Stabilization parameter 1
        ##  STABILIZATION functional terms
        # mass conservation residual
        mass_residual  = s_gauss # Source Term
        mass_residual += -Decay_gauss*f_gauss # Decay_gauss factor
        # Filter subscales
        filter_subscale = tau1*mass_residual
        # ASGS stabilization
        # some of the terms in rv_stab are taken with sign (+) since they derive from an integration by parts 
        rv_stab = -Decay_gauss*q_gauss.transpose()*filter_subscale # stab decay residual
        rv += rv_stab # save STABILIZATION residual into TOTAL element residual
    # END HANDLE STABILIZATION

    ## Define DOFs and test function vectors
    dofs = sympy.zeros(nnodes, 1)
    testfunc = sympy.zeros(nnodes, 1)
    for i in range(nnodes):
        # Filter DOFs and test functions
        dofs[i] = f[i,0]
        testfunc[i] = q[i,0]

    # COMPUTE RHS AND LHS
    print(f"Computing {dim}D{nnodes}N RHS Gauss point contribution\n")
    rhs = Compute_RHS(rv.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(gauss_weight*rhs, "rRHS", mode, assignment_op='+=')
    print(f"Computing {dim}D{nnodes}N LHS Gauss point contribution\n")
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications) # Compute the LHS (considering stress as C*(B*v) to derive w.r.t. v)
    lhs_out = OutputMatrix_CollectingFactors(gauss_weight*lhs, "rLHS", mode, assignment_op='+=')
    # print(OutputMatrix_CollectingFactors(lhs,"lhs",mode))
    ## Replace the computed RHS and LHS in the template outstring
    outstring = outstring.replace(f"//substitute_lhs_{dim}D{nnodes}N", lhs_out)
    outstring = outstring.replace(f"//substitute_rhs_{dim}D{nnodes}N", rhs_out)
    ## | END ELEMENT DEFINITION |
    ## |---------------------------------------------|
    

## Write the modified template
print("Writing output file \'" + output_filename + "\'")
out = open(output_filename,'w')
out.write(outstring)
out.close()

