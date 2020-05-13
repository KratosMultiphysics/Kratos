from sympy import *
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

## Settings explanation
# DIMENSION TO COMPUTE:
# This symbolic generator is valid for both 2D and 3D cases. Since the element has been programed with a dimension template in Kratos,
# it is advised to set the dim_to_compute flag as "Both". In this case the generated .cpp file will contain both 2D and 3D implementations.

## Symbolic generation settings
do_simplifications = False
dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
ASGS_stabilization = False          # Consider ASGS stabilization terms
mode = "c"                          # Output mode to a c++ file

if (dim_to_compute == "2D"):
    dim_vector = [2]
elif (dim_to_compute == "3D"):
    dim_vector = [3]
elif (dim_to_compute == "Both"):
    dim_vector = [2,3]

## Read the template file and set the output filename accordingly
template_filename = "symbolic_eulerian_convection_diffusion_explicit_cpp_template.cpp"

if template_filename == "symbolic_eulerian_convection_diffusion_explicit_cpp_template.cpp":
    output_filename = "symbolic_eulerian_convection_diffusion_explicit.cpp"
else:
    err_msg = "Wrong template_filename provided. Must be (template --> output):\n"
    err_msg +=  "\t- symbolic_eulerian_convection_diffusion_explicit_cpp_template.cpp --> symbolic_eulerian_convection_diffusion_explicit.cpp"
    raise Exception(err_msg)

## Initialize the outstring to be filled with the template .cpp file
templatefile = open(template_filename)
outstring = templatefile.read()

for dim in dim_vector:
    if(dim == 2):
        nnodes = 3
    elif(dim == 3):
        nnodes = 4

    impose_partion_of_unity = False
    N,DN = DefineShapeFunctions(nnodes, dim, impose_partion_of_unity)

    ## Unknown fields definition
    phi = DefineVector('phi',nnodes) # scalar unknown

    ## Test functions definition
    w = DefineMatrix('w',nnodes,dim) # vector unknown field test function (not needed)
    q = DefineVector('q',nnodes)     # scalar unknown field test function

    ## Other data definitions
    f = DefineVector('f',nnodes)     # forcing term

    ## Other simbols definition
    k = Symbol('k',positive= True)   # diffusion coefficient
    v = DefineMatrix('v',nnodes,dim) # convective velocity

    ## Data interpolation to the Gauss points
    f_gauss = f.transpose()*N
    w_gauss = w.transpose()*N
    q_gauss = q.transpose()*N
    v_gauss = v.transpose()*N
    phi_gauss = phi.transpose()*N

    ## Gradient and divergence computation
    grad_w = DfjDxi(DN,w)
    grad_q = DfjDxi(DN,q)
    grad_phi = DfjDxi(DN,phi)
    div_w = div(DN,w)

    ## Compute galerkin functional
    # Diffusion functional
    rhs_forcing = q_gauss.transpose() * f_gauss
    rhs_diffusion = - k * grad_phi.transpose() * grad_q
    rhs_convective = - q_gauss * (v_gauss.transpose() * grad_phi)
    # rhs_convective = - phi_gauss * (v_gauss.transpose() * grad_q) # gives equivalent result
    rhs_galerkin = rhs_forcing + rhs_diffusion + rhs_convective

    ##  Stabilization functional terms

    # Mass conservation residual

    # Compute the ASGS stabilization terms using the momentum and mass conservation residuals above

    ## Add the stabilization terms to the original residual terms
    if (ASGS_stabilization):
        res = rhs_galerkin + rhs_stab
    else:
        res = rhs_galerkin

    ## Define DOFs and test function vectors
    dofs = Matrix( zeros(nnodes*(1), 1) ) # 1 because scalar unknown
    testfunc = Matrix( zeros(nnodes*(1), 1) ) # 1 because scalar unknown

    for i in range(0,nnodes):
        # phi DOFs and test functions
        dofs[i*(1)] = phi[i,0]
        testfunc[i*(1)] = q[i,0]

    ## Compute LHS and RHS
    rhs = Compute_RHS(res.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)

    # Compute LHS (RHS(residual) differenctiation w.r.t. the DOFs)
    lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications)
    lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)

    ## Replace the computed RHS and LHS in the template outstring
    if(dim == 2):
        outstring = outstring.replace("//substitute_lhs_2D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)
    elif(dim == 3):
        outstring = outstring.replace("//substitute_lhs_3D", lhs_out)
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

## Write the modified template
out = open(output_filename,'w')
out.write(outstring)
out.close()