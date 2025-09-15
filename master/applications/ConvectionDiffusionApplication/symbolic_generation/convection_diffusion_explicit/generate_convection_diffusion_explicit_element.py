from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *
from sympy import *

def DefineShapeFunctionsMatrix(dim, n_nodes, n_gauss):
    mat_N = DefineMatrix('mat_N', n_gauss, n_nodes)
    if dim == 2:
        if n_gauss == 1:
            mat_N[0,0] = 1.0 / 3.0
            mat_N[0,1] = 1.0 / 3.0
            mat_N[0,2] = 1.0 / 3.0
        elif n_gauss == 3:
            mat_N[0,0] = 2.0 / 3.0
            mat_N[0,1] = 1.0 / 6.0
            mat_N[0,2] = 1.0 / 6.0
            mat_N[1,0] = 1.0 / 6.0
            mat_N[1,1] = 2.0 / 3.0
            mat_N[1,2] = 1.0 / 6.0
            mat_N[2,0] = 1.0 / 6.0
            mat_N[2,1] = 1.0 / 6.0
            mat_N[2,2] = 2.0 / 3.0
        else:
            err_msg = "Invalid quadrature for dimension " + str(dim) + " and number of Gauss points " + str(n_gauss) + "."
    elif dim == 3:
        if n_gauss == 1:
            mat_N[0,0] = 1.0 / 4.0
            mat_N[0,1] = 1.0 / 4.0
            mat_N[0,2] = 1.0 / 4.0
            mat_N[0,3] = 1.0 / 4.0
        elif n_gauss == 4:
            mat_N[0,0] = 0.58541020
            mat_N[0,1] = 0.13819660
            mat_N[0,2] = 0.13819660
            mat_N[0,3] = 0.13819660
            mat_N[1,0] = 0.13819660
            mat_N[1,1] = 0.58541020
            mat_N[1,2] = 0.13819660
            mat_N[1,3] = 0.13819660
            mat_N[2,0] = 0.13819660
            mat_N[2,1] = 0.13819660
            mat_N[2,2] = 0.58541020
            mat_N[2,3] = 0.13819660
            mat_N[3,0] = 0.13819660
            mat_N[3,1] = 0.13819660
            mat_N[3,2] = 0.13819660
            mat_N[3,3] = 0.58541020
        else:
            err_msg = "Invalid quadrature for dimension " + str(dim) + " and number of Gauss points " + str(n_gauss) + "."
    else:
        err_msg = "Invalid dimension " + str(dim) + "."
        raise Exception(err_msg)

    return mat_N

# Symbolic generation settings
do_simplifications = False
dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
stabilization = True                # Consider ASGS or OSS stabilization. By default simple ASGS is used
OSS_stabilization = True            # Requires stabilization to be true
dynamic_subscales = False            # Consider subscale dynamic
mode = "c"                          # Output mode to a c++ file

if (dim_to_compute == "2D"):
    dim_vector = [2]
elif (dim_to_compute == "3D"):
    dim_vector = [3]
elif (dim_to_compute == "Both"):
    dim_vector = [2,3]

# Initialize the outstring to be filled with the template .cpp file
if dynamic_subscales:
    template_filename = "d_convection_diffusion_explicit_cpp_template.cpp"
else:
    template_filename = "qs_convection_diffusion_explicit_cpp_template.cpp"
templatefile = open(template_filename)
outstring = templatefile.read()

# Set the output filename
if template_filename == "d_convection_diffusion_explicit_cpp_template.cpp":
    output_filename = "d_convection_diffusion_explicit.cpp"
elif template_filename == "qs_convection_diffusion_explicit_cpp_template.cpp":
    output_filename = "qs_convection_diffusion_explicit.cpp"
else:
    err_msg = "Wrong template_filename provided. Must be (template --> output):\n"
    err_msg +=  "\t- d_convection_diffusion_explicit_cpp_template.cpp --> d_convection_diffusion_explicit.cpp\n"
    err_msg +=  "\t- qs_convection_diffusion_explicit_cpp_template.cpp --> qs_convection_diffusion_explicit.cpp"
    raise Exception(err_msg)

for dim in dim_vector:
    # Shape functions and Gauss points settings
    if(dim == 2):
        nnodes = 3
        ngauss = 3
    elif(dim == 3):
        nnodes = 4
        ngauss = 4

    DN = DefineMatrix('DN', nnodes, dim)
    mat_N = DefineShapeFunctionsMatrix(dim, nnodes, ngauss)

    # Unknown fields definition
    phi = DefineVector('phi',nnodes)                  # scalar unknown
    phi_old = DefineVector('phi_old',nnodes)          # scalar unknown previous time step
    phi_subscale_gauss = DefineVector('phi_subscale_gauss',ngauss) # scalar unknown on subscale space, defined on gauss integration points
    phi_acceleration_old = DefineVector('phi_acceleration_old',nnodes)          # acceleration scalar unknown previous time step

    # Test functions definition
    w = DefineMatrix('w',nnodes,dim)   # vector unknown field test function (not needed)
    q = DefineVector('q',nnodes)       # scalar unknown field test function

    # Other data definitions
    f = DefineVector('f',nnodes)                                    # forcing term
    alpha = Symbol('alpha',positive= True)                          # diffusion coefficient
    v = DefineMatrix('v',nnodes,dim)                                # convective velocity
    tau = DefineVector('tau',ngauss)                                # stabilization coefficient
    delta_time = Symbol('delta_time',positive= True)                # time current time step
    explicit_step_coefficient = Symbol('explicit_step_coefficient') # coefficient for explicit scheme
    prj = DefineVector('prj',nnodes)                                # OSS projection term

    # Loop and accumulate the residual in each Gauss point
    res_tot = Matrix(zeros(1,1))
    res_OSS_tot = Matrix(zeros(1,1))

    for i_gauss in range(ngauss):
        print("\tGauss point: " + str(i_gauss))

        # Get Gauss point geometry data
        N = DefineVector('N', nnodes)
        for i in range(nnodes):
            N[i] = mat_N[i_gauss,i]

        ########## Estimate variables on gauss point ##########

        # Data interpolation to the Gauss points
        f_gauss = f.transpose()*N
        w_gauss = w.transpose()*N
        q_gauss = q.transpose()*N
        v_gauss = v.transpose()*N
        phi_gauss = phi.transpose()*N
        phi_old_gauss = phi_old.transpose()*N
        phi_acceleration_old_gauss = phi_acceleration_old.transpose()*N
        prj_gauss = prj.transpose()*N

        # Gradients and divergences computation
        grad_w = DfjDxi(DN,w)
        grad_q = DfjDxi(DN,q)
        grad_phi = DfjDxi(DN,phi)
        grad_f = DfjDxi(DN,f)
        div_w = ones(1,1)*sum([DN[i]*w[i] for i in range (0,len(DN))])
        div_v = ones(1,1)*sum([DN[i]*v[i] for i in range (0,len(DN))])

        ########## Estimate residuals on gauss point ##########

        ##### Galerkin functional terms #####
        rhs_forcing = q_gauss.transpose() * f_gauss
        rhs_diffusion = - alpha * grad_phi.transpose() * grad_q
        rhs_convective_1 = - q_gauss * (v_gauss.transpose() * grad_phi)
        rhs_convective_2 = - q_gauss * phi_gauss * div_v
        rhs_galerkin = rhs_forcing + rhs_diffusion + rhs_convective_1 + rhs_convective_2

        #####  Stabilization ASGS/OSS functional terms #####
        # ASGS/OSS Convective term
        rhs_stab_1_forcing = tau[i_gauss] * (v_gauss.transpose() * grad_q) * f_gauss
        rhs_stab_1_mass = - tau[i_gauss] * (grad_q.transpose() * v_gauss) * N.transpose() * (phi-phi_old)*explicit_step_coefficient
        rhs_stab_1_convection_1 = - tau[i_gauss] * (v_gauss.transpose() * grad_q) * (v_gauss.transpose() * grad_phi)
        rhs_stab_1_convection_2 = - tau[i_gauss] * (v_gauss.transpose() * grad_q) * phi_gauss * div_v
        rhs_stab_1 = rhs_stab_1_forcing + rhs_stab_1_convection_1 + rhs_stab_1_convection_2 + rhs_stab_1_mass
        # OSS term of convective term
        rhs_stab_1_oss = tau[i_gauss] * (v_gauss.transpose() * grad_q) * prj_gauss
        if OSS_stabilization:
            rhs_stab_1 += rhs_stab_1_oss
        # Dynamic term of convective term
        rhs_stab_1_dynamic = tau[i_gauss] * (v_gauss.transpose() * grad_q) * phi_subscale_gauss[i_gauss]/delta_time
        if dynamic_subscales:
            rhs_stab_1 += rhs_stab_1_dynamic

        # ASGS/OSS dynamic term
        rhs_stab_2_forcing = - tau[i_gauss] * q_gauss.transpose() * f_gauss
        rhs_stab_2_mass = tau[i_gauss] * q_gauss.transpose() * N.transpose() * (phi-phi_old)*explicit_step_coefficient
        rhs_stab_2_convection_1 = tau[i_gauss] * q_gauss * (v_gauss.transpose() * grad_phi)
        rhs_stab_2_convection_2 = tau[i_gauss] * q_gauss * phi_gauss * div_v
        rhs_stab_2_diffusion = tau[i_gauss] * alpha * grad_phi.transpose() * grad_q
        rhs_stab_2_subgrid_old = - tau[i_gauss] * q_gauss.transpose() * phi_subscale_gauss[i_gauss]
        rhs_stab_2_mass_subgrid_old = q_gauss.transpose() * phi_subscale_gauss[i_gauss] / delta_time
        rhs_stab_2 = rhs_stab_2_forcing + rhs_stab_2_mass + rhs_stab_2_convection_1 + rhs_stab_2_convection_2 + rhs_stab_2_diffusion + rhs_stab_2_subgrid_old + rhs_stab_2_mass_subgrid_old
        # OSS term of mass term
        rhs_stab_2_oss = - tau[i_gauss] * q_gauss.transpose() * prj_gauss
        if OSS_stabilization:
            rhs_stab_2 += rhs_stab_2_oss

        # Compute ASGS/OSS stabilization functional
        rhs_stabilization = rhs_stab_1
        # Eventually add ASGS/OSS dynamic term
        if dynamic_subscales:
            rhs_stabilization = rhs_stabilization + rhs_stab_2

        ##### OSS step #####
        # with lhs we refer to the fact we take the strong equation on the left side
        lhs_OSS_forcing = - q_gauss.transpose() * f_gauss
        lhs_OSS_mass = q_gauss.transpose() * (N.transpose() * (phi-phi_old)*explicit_step_coefficient)
        lhs_OSS_mass_subscale = - q_gauss.transpose() * (phi_subscale_gauss[i_gauss]/delta_time)
        lhs_OSS_diffusion = alpha * grad_phi.transpose() * grad_q
        lhs_OSS_convective_1 = q_gauss * (v_gauss.transpose() * grad_phi)
        lhs_OSS_convective_2 = q_gauss * phi_gauss * div_v
        res_OSS = lhs_OSS_forcing + lhs_OSS_mass + lhs_OSS_diffusion + lhs_OSS_convective_1 + lhs_OSS_convective_2
        if dynamic_subscales:
            res_OSS = res_OSS + lhs_OSS_mass_subscale

        # Add the stabilization terms to the original residual terms
        if (stabilization):
            res = rhs_galerkin + rhs_stabilization
        else:
            res = rhs_galerkin

        # Accumulate in the total value
        res_tot += res
        res_OSS_tot += res_OSS

    # Define DOFs and test function matrices
    dofs = Matrix( zeros(nnodes*(1), 1) ) # 1 because scalar unknown
    testfunc = Matrix( zeros(nnodes*(1), 1) ) # 1 because scalar unknown
    for i in range(0,nnodes):
        # phi DOFs and test functions
        dofs[i*(1)] = phi[i,0]
        testfunc[i*(1)] = q[i,0]

    # Compute RHS
    rhs = Compute_RHS(res_tot.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)

    # Replace the computed RHS and LHS in the template outstring
    if(dim == 2):
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)
    elif(dim == 3):
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

    # Compute OSS functional
    if OSS_stabilization:
        rhs = Compute_RHS(res_OSS_tot.copy(), testfunc, do_simplifications)
        rhs_out = OutputVector_CollectingFactors(rhs, "rhs", mode)

        # Replace the computed OSS in the template outstring
        if(dim == 2):
            outstring = outstring.replace("//substitute_oss_2D", rhs_out)
        elif(dim == 3):
            outstring = outstring.replace("//substitute_oss_3D", rhs_out)

    # Substitute the shape function gradients container accesses
    for i_node in range(nnodes):
        for j_dim in range(dim):
            to_substitute = 'DN(' + str(i_node) + ',' + str(j_dim) + ')'
            substituted_value = 'DN_DX_' + str(i_node) + '_' + str(j_dim)
            outstring = outstring.replace(to_substitute, substituted_value)

# Write the modified template
out = open(output_filename,'w')
out.write(outstring)
out.close()
