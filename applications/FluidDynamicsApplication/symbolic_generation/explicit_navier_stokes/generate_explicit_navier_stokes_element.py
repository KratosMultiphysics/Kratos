from sympy import *
from KratosMultiphysics import *
from KratosMultiphysics.sympy_fe_utilities import *

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
mode = "c"
do_simplifications = False
dim_to_compute = "Both"             # Spatial dimensions to compute. Options:  "2D","3D","Both"
ASGS_stabilization = False          # Consider ASGS stabilization terms
OSS_stabilization = False           # Consider OSS stabilization terms
formulation = "NavierStokes"        # Element type. Options: "NavierStokes"

if formulation == "NavierStokes":
    convective_term = True
    linearisation = "Picard" # Convective term linearisation type. Options: "Picard", "FullNR"
    output_filename = "qs_navier_stokes_explicit.cpp"
    template_filename = "symbolic_explicit_qs_navier_stokes_cpp_template.cpp"
else:
    err_msg = "Wrong formulation. Given \'" + formulation + "\'. Available option is \'NavierStokes\'."
    raise Exception(err_msg)

info_msg = "\n"
info_msg += "Element generator settings:\n"
info_msg += "\t - Element type: " + formulation + "\n"
info_msg += "\t - Dimension: " + dim_to_compute + "\n"
info_msg += "\t - ASGS stabilization: " + str(ASGS_stabilization) + "\n"
info_msg += "\t - OSS stabilization: " + str(OSS_stabilization) + "\n"
print(info_msg)

if formulation == "NavierStokes":
    if (dim_to_compute == "2D"):
        dim_vector = [2]
        n_nodes_vector = [3] # tria
    elif (dim_to_compute == "3D"):
        dim_vector = [3]
        n_nodes_vector = [4] # tet
    elif (dim_to_compute == "Both"):
        dim_vector = [2, 3] # tria, tet

# Initialize the outstring to be filled with the template .cpp file
print("Reading template file \'"+ template_filename + "\'\n")
templatefile = open(template_filename)
outstring = templatefile.read()

for dim in dim_vector:

    # Shape functions and Gauss points settings
    if(dim == 2):
        n_nodes = 3
        n_gauss = 3
    elif(dim == 3):
        n_nodes = 4
        n_gauss = 4

    DN = DefineMatrix('DN', n_nodes, dim)
    mat_N = DefineShapeFunctionsMatrix(dim, n_nodes, n_gauss)

    # Unknown fields definition
    v = DefineMatrix('v',n_nodes,dim)            # Current step velocity (v(i,j) refers to velocity of node i component j)
    vn = DefineMatrix('vn',n_nodes,dim)          # Previous step velocity
    p = DefineVector('p',n_nodes)                # Pressure
    pn = DefineVector('pn',n_nodes)              # Previous step pressure

    # Test functions definition
    w = DefineMatrix('w',n_nodes,dim)            # Velocity field test function
    q = DefineVector('q',n_nodes)                # Pressure field test function

    # Other data definitions
    f = DefineMatrix('f',n_nodes,dim)            # Forcing term
    dt  = Symbol('dt', positive = True)          # Time increment
    rho = Symbol('rho', positive = True)         # Density
    nu  = Symbol('nu', positive = True)          # Kinematic viscosity (mu/rho)
    mu  = Symbol('mu', positive = True)          # Dynamic viscosity
    h = Symbol('h', positive = True)
    dyn_tau = Symbol('dyn_tau', positive = True)
    stab_c1 = Symbol('stab_c1', positive = True)
    stab_c2 = Symbol('stab_c2', positive = True)

    res_tot = Matrix(zeros(1,1))

    for i_gauss in range(n_gauss):
        print("\tGauss point: " + str(i_gauss))

        # Get Gauss point geometry data
        N = DefineVector('N', n_nodes)
        for i in range(n_nodes):
            N[i] = mat_N[i_gauss,i]

        # Data interpolation to the Gauss points
        f_gauss = f.transpose()*N
        v_gauss = v.transpose()*N
        p_gauss = p.transpose()*N
        w_gauss = w.transpose()*N
        q_gauss = q.transpose()*N

        # Gradients computation (fluid dynamics gradient)
        grad_w = DfjDxi(DN,w)
        grad_q = DfjDxi(DN,q)
        grad_p = DfjDxi(DN,p)
        grad_v = DfjDxi(DN,v)
        div_v = ones(1,1)*sum([DN[i]*v[i] for i in range (0,len(DN))])
        div_w = ones(1,1)*sum([DN[i]*w[i] for i in range (0,len(DN))])
        # TODO: check following is correct
        grad_sym_v = grad_sym_voigtform(DN,v) # Symmetric gradient of v in Voigt notation
        grad_sym_w = grad_sym_voigtform(DN,w) # Symmetric gradient of w in Voigt notation
        # Recall that the grad(w):grad(v) contraction equals grad_sym(w)*grad_sym(v) in Voigt notation since they are symmetric tensors.

        # Convective velocity definition
        if convective_term:
            if (linearisation == "Picard"):
                vconv = DefineMatrix('vconv',n_nodes,dim) # Convective velocity defined a symbol
            elif (linearisation == "FullNR"):
                vmesh = DefineMatrix('vmesh',n_nodes,dim) # Mesh velocity
                vconv = v - vmesh                         # Convective velocity defined as a velocity dependent variable
            else:
                raise Exception("Wrong linearisation \'" + linearisation + "\' selected. Available options are \'Picard\' and \'FullNR\'.")
            vconv_gauss = vconv.transpose()*N
        if convective_term:
            div_vconv = div(DN,vconv)
            div_vconv = ones(1,1)*sum([DN[i]*vconv[i] for i in range (0,len(DN))])
        if convective_term:
            convective_term_gauss = (vconv_gauss.transpose()*grad_v)

        # Compute the stabilization parameters
        if convective_term:
            stab_norm_a = 0.0
            for i in range(0, dim):
                stab_norm_a += vconv_gauss[i]**2
            stab_norm_a = sqrt(stab_norm_a)
            tau1 = 1.0/((rho*dyn_tau)/dt + (stab_c2*rho*stab_norm_a)/h + (stab_c1*mu)/(h*h)) # Stabilization parameter 1
            tau2 = mu + (stab_c2*rho*stab_norm_a*h)/stab_c1                                  # Stabilization parameter 2
        else:
            tau1 = 1.0/((rho*dyn_tau)/dt + (stab_c1*mu)/(h*h)) # Stabilization parameter 1
            tau2 = (h*h) / (stab_c1 * tau1)                    # Stabilization parameter 2

        # Compute functionals

        # Compute galerkin functional
        # Navier-Stokes functional
        res_galerkin = w_gauss.transpose()*f_gauss - nu*grad_sym_w.transpose()*grad_sym_v + div_w*p_gauss - q_gauss*div_v
        if convective_term:
            res_galerkin -= w_gauss.transpose()*convective_term_gauss.transpose()

        # Add the stabilization terms to the original residual terms
        if (ASGS_stabilization):
            res = res_galerkin + res_stabilization
        else:
            res = res_galerkin

        # Accumulate in the total residual
        res_tot += res

    # Define DOFs and test function vectors
    dofs = Matrix( zeros(n_nodes*(dim+1), 1) )
    testfunc = Matrix( zeros(n_nodes*(dim+1), 1) )

    for i in range(0,n_nodes):

        # Velocity DOFs and test functions
        for k in range(0,dim):
            dofs[i*(dim+1)+k] = v[i,k]
            testfunc[i*(dim+1)+k] = w[i,k]

        # Pressure DOFs and test functions
        dofs[i*(dim+1)+dim] = p[i,0]
        testfunc[i*(dim+1)+dim] = q[i,0]

    # Compute LHS and RHS
    print("Computing " + str(dim) + "D RHS Gauss point contribution\n")
    rhs = Compute_RHS(res_tot.copy(), testfunc, do_simplifications)
    rhs_out = OutputVector_CollectingFactors(rhs, "rRightHandSideBoundedVector", mode)

    # print("Computing " + str(dim) + "D LHS Gauss point contribution\n")
    # lhs = Compute_LHS(rhs, testfunc, dofs, do_simplifications)
    # lhs_out = OutputMatrix_CollectingFactors(lhs, "lhs", mode)

    # Replace the computed RHS and LHS in the template outstring
    if(dim == 2):
        outstring = outstring.replace("//substitute_rhs_2D", rhs_out)
    elif(dim == 3):
        outstring = outstring.replace("//substitute_rhs_3D", rhs_out)

    # Substitute the shape function gradients container accesses
    for i_node in range(n_nodes):
        for j_dim in range(dim):
            to_substitute = 'DN(' + str(i_node) + ',' + str(j_dim) + ')'
            substituted_value = 'DN_DX_' + str(i_node) + '_' + str(j_dim)
            outstring = outstring.replace(to_substitute, substituted_value)

# Write the modified template
print("Writing output file \'" + output_filename + "\'")
out = open(output_filename,'w')
out.write(outstring)
out.close()
