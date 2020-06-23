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
ASGS_stabilization = True           # Consider ASGS stabilization terms
dynamic_subscales = True            # Consider subscale dynamic
OSS_stabilization = True            # Consider OSS stabilization terms (together with ASGS: requires ASGS stabilization to be true)
mode = "c"                          # Output mode to a c++ file

if (dim_to_compute == "2D"):
    dim_vector = [2]
elif (dim_to_compute == "3D"):
    dim_vector = [3]
elif (dim_to_compute == "Both"):
    dim_vector = [2,3]

## Read the template file and set the output filename accordingly
if dynamic_subscales:
    template_filename = "symbolic_dynamic_eulerian_convection_diffusion_explicit_cpp_template.cpp"
else:
    template_filename = "symbolic_quasi_static_eulerian_convection_diffusion_explicit_cpp_template.cpp"

if template_filename == "symbolic_dynamic_eulerian_convection_diffusion_explicit_cpp_template.cpp":
    output_filename = "symbolic_dynamic_eulerian_convection_diffusion_explicit.cpp"
elif template_filename == "symbolic_quasi_static_eulerian_convection_diffusion_explicit_cpp_template.cpp":
    output_filename = "symbolic_quasi_static_eulerian_convection_diffusion_explicit.cpp"
else:
    err_msg = "Wrong template_filename provided. Must be (template --> output):\n"
    err_msg +=  "\t- symbolic_dynamic_eulerian_convection_diffusion_explicit_cpp_template.cpp --> symbolic_dynamic_eulerian_convection_diffusion_explicit.cpp\n"
    err_msg +=  "\t- symbolic_quasi_static_eulerian_convection_diffusion_explicit_cpp_template.cpp --> symbolic_quasi_static_eulerian_convection_diffusion_explicit.cpp"
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
    phi = DefineVector('phi',nnodes)                  # scalar unknown
    phi_old = DefineVector('phi_old',nnodes)          # scalar unknown previous time step
    phi_subscale_gauss = Symbol('phi_subscale_gauss') # scalar unknown on subscale space, defined on gauss integration points

    ## Test functions definition
    w = DefineMatrix('w',nnodes,dim)   # vector unknown field test function (not needed)
    q = DefineVector('q',nnodes)       # scalar unknown field test function

    ## Other data definitions
    f = DefineVector('f',nnodes)       # forcing term
    k = Symbol('k',positive= True)     # diffusion coefficient
    v = DefineMatrix('v',nnodes,dim)   # convective velocity
    tau = Symbol('tau',positive= True) # stabilization coefficient
    delta_time = Symbol('delta_time',positive= True) # time current time step
    RK_time_coefficient = Symbol('RK_time_coefficient',positive= True) # time coefficient for RK scheme
    prj = DefineVector('prj',nnodes)       # OSS projection term

    ## Data interpolation to the Gauss points
    f_gauss = f.transpose()*N
    w_gauss = w.transpose()*N
    q_gauss = q.transpose()*N
    v_gauss = v.transpose()*N
    phi_gauss = phi.transpose()*N
    phi_old_gauss = phi_old.transpose()*N
    prj_gauss = prj.transpose()*N

    ## Gradient and divergence computation
    grad_w = DfjDxi(DN,w)
    grad_q = DfjDxi(DN,q)
    grad_phi = DfjDxi(DN,phi)
    div_w = div(DN,w)
    div_v = div(DN,v)
    grad_f = DfjDxi(DN,f)

    ## Galerkin functional terms
    rhs_forcing = q_gauss.transpose() * f_gauss
    rhs_diffusion = - k * grad_phi.transpose() * grad_q
    rhs_convective_1 = - q_gauss * (v_gauss.transpose() * grad_phi)
    rhs_convective_2 = - q_gauss * phi_gauss * div_v
    rhs_galerkin = rhs_forcing + rhs_diffusion + rhs_convective_1 + rhs_convective_2

    ##  Stabilization ASGS functional terms
    # Convective term
    rhs_stab_1_forcing = tau * (v_gauss.transpose() * grad_q) * f_gauss
    rhs_stab_1_mass = - tau * (grad_q.transpose() * v_gauss) * N.transpose() * (phi-phi_old)/(RK_time_coefficient*delta_time)
    rhs_stab_1_convection_1 = - tau * (v_gauss.transpose() * grad_q) * (v_gauss.transpose() * grad_phi)
    rhs_stab_1_convection_2 = - tau * (v_gauss.transpose() * grad_q) * phi_gauss * div_v
    rhs_stab_1 = rhs_stab_1_forcing + rhs_stab_1_convection_1 + rhs_stab_1_convection_2 + rhs_stab_1_mass
    # OSS term of convective term
    rhs_stab_1_oss = tau * (v_gauss.transpose() * grad_q) * prj_gauss
    if OSS_stabilization:
        rhs_stab_1 += rhs_stab_1_oss
    # dynamic term of convective term
    rhs_stab_1_dynamic = tau * (v_gauss.transpose() * grad_q) * phi_subscale_gauss/delta_time
    if dynamic_subscales:
        rhs_stab_1 += rhs_stab_1_dynamic

    # Mass conservation residual
    # mass term
    rhs_stab_2_forcing = - q_gauss.transpose() * f_gauss
    rhs_stab_2_mass = q_gauss.transpose() * N.transpose() * (phi-phi_old)/(RK_time_coefficient*delta_time)
    rhs_stab_2_convection_1 = q_gauss * (v_gauss.transpose() * grad_phi)
    rhs_stab_2_convection_2 = q_gauss * phi_gauss * div_v
    rhs_stab_2_diffusion = k * grad_phi.transpose() * grad_q
    rhs_stab_2_mass_subgrid_old = (1/tau) * q_gauss.transpose() * phi_subscale_gauss
    rhs_stab_2 = rhs_stab_2_forcing + rhs_stab_2_mass + rhs_stab_2_convection_1 + rhs_stab_2_convection_2 + rhs_stab_2_diffusion + rhs_stab_2_mass_subgrid_old
    # OSS term of mass term
    rhs_stab_2_oss = - q_gauss.transpose() * prj_gauss
    if OSS_stabilization:
        rhs_stab_2 += rhs_stab_2_oss

    # Compute the ASGS stabilization terms using the momentum and mass conservation residuals above
    rhs_stabilization = rhs_stab_1
    if dynamic_subscales:
        rhs_stabilization += rhs_stab_2

    ## Stabilization OSS funtional terms
    # with lhs we refer to the fact we take the strong equation on the left side
    lhs_OSS_forcing = - q_gauss.transpose() * f_gauss
    lhs_OSS_mass = q_gauss.transpose() * (N.transpose() * (phi-phi_old)/(RK_time_coefficient*delta_time))
    lhs_OSS_mass_subscale = - q_gauss.transpose() * (phi_subscale_gauss/delta_time)
    lhs_OSS_diffusion = k * grad_phi.transpose() * grad_q
    lhs_OSS_convective_1 = q_gauss * (v_gauss.transpose() * grad_phi)
    lhs_OSS_convective_2 = q_gauss * phi_gauss * div_v
    res_OSS = lhs_OSS_forcing + lhs_OSS_mass + lhs_OSS_diffusion + lhs_OSS_convective_1 + lhs_OSS_convective_2
    if dynamic_subscales:
        res_OSS = res_OSS + lhs_OSS_mass_subscale

    ## Add the stabilization terms to the original residual terms
    if (ASGS_stabilization):
        res = rhs_galerkin + rhs_stabilization
    else:
        res = rhs_galerkin

    ## Define DOFs and test function vectors
    dofs = Matrix( zeros(nnodes*(1), 1) ) # 1 because scalar unknown
    testfunc = Matrix( zeros(nnodes*(1), 1) ) # 1 because scalar unknown

    for i in range(0,nnodes):
        # phi DOFs and test functions
        dofs[i*(1)] = phi[i,0]
        testfunc[i*(1)] = q[i,0]

    ## Compute RHS
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

    if OSS_stabilization:
        # Compute OSS residual
        rhs = Compute_RHS(res_OSS.copy(), testfunc, do_simplifications)
        rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)

        ## Replace the computed OSS in the template outstring
        if(dim == 2):
            outstring = outstring.replace("//substitute_oss_2D", rhs_out)
        elif(dim == 3):
            outstring = outstring.replace("//substitute_oss_3D", rhs_out)

## Write the modified template
out = open(output_filename,'w')
out.write(outstring)
out.close()
