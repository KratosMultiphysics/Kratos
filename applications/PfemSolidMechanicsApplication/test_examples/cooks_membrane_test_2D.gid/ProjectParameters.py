domain_size = 2

#Problem Data
#####################################

ProblemType = "Mechanical"
NumberofThreads = 4
Solution_method = "Newton-Raphson"                
SolverType  = "QuasiStaticSolver"
time_step = 1
nsteps = 1

#Solver Data
#####################################

class SolverSettings:
    solver_type  = "mechanical_solver"
    domain_size  = 2
    RotationDofs = False
    PressureDofs = True
    ReformDofSetAtEachStep = True
    compute_reactions = True
    line_search = "False"
    scheme_type = "QuasiStaticSolver"
    convergence_criterion = "Residual_Criteria"
    displacement_relative_tolerance = 1.00000e-09
    displacement_absolute_tolerance = 1.00000e-09
    residual_relative_tolerance = 1.00000e-09
    residual_absolute_tolerance = 1.00000e-09
    max_iteration =   30

    class linear_solver_config:
        solver_type   = "Super LU"
        scaling       = False
        tolerance     = 1e-7
        max_iteration = 5000 #300 
        verbosity     = 0
        is_symmetric  = False  
        
        #Pastix Iterative Solver:
        gmres_krylov_space_dimension = 100
        ilu_level_of_fill            = 3 #5

        #GMRES or CG:
        preconditioner_type          = "None"

        #Deflated CG:
        assume_constant_structure    = True
        max_reduced_size             = 1000

        #AMG: (requires ResidualBasedBlockBuilderAndSolver )
        smoother_type  = "ILU0" #"DAMPED_JACOBI"
        krylov_type    = "GMRES"
        

#Meshing Data
#####################################

RemeshDomains  = "False"
MeshingElement = "SpatialLagrangianUPElement2D3N"
RefineDomains  = "False"
MeshConditions = []
Conditions1 = {"Subdomain":1, "Remesh":1, "Constrained":1, "Refine":1, "MeshSmoothing":1, "JacobiSmoothing":1, "MeshElement":"SpatialLagrangianUPElement2D3N"}
MeshConditions.append(Conditions1)

BoxCenter = []
BoxCenter.append(0)
BoxCenter.append(0)
BoxCenter.append(0)

BoxVelocity = [] 
BoxVelocity.append(0)
BoxVelocity.append(0)
BoxVelocity.append(0)

#set mesh modeler configuration
class mesh_modeler_config:
    number_domains = 1
    size_scale = 1
    critical_radius = 0.025
    critical_dissipation = 100
    reference_error = 2
    offset_factor = 0.0001                             
    mesh_conditions = MeshConditions
    box_refinement_only = "False"
    box_center = BoxCenter
    box_velocity = BoxVelocity
    box_radius = 0
    remesh_frequency = 1

#Contact Data
#####################################

FindContacts = "False"

class contact_modeler_config:
    contact_condition        = "ContactDomainLM2DCondition"
    constrained_contact      = "True"
    friction_active          = "False"
    mu_static                = 0.3
    mu_dynamic               = 0.2
    offset_factor            = 0.0001
    penalty_parameter        = 1000
    stability_parameter      = 0.01
    contact_search_frequency = 1


ToolVelocity = []
ToolVelocity.append(0)
ToolVelocity.append(0)
ToolVelocity.append(0)
        
ToolCenter = []
ToolCenter.append(0)
ToolCenter.append(0)
ToolCenter.append(0)


#set rigid wall configuration
class rigid_wall_config:
    rigid_wall         = "False"
    size_scale         = 1
    tip_radius         = 0.0004
    rake_angle         = 5
    clearance_angle    = 5
    young_modulus      = 210e9
    penalty_parameter  = 2.1
    center             = ToolCenter
    velocity           = ToolVelocity


#Constraints Data
#####################################

Incremental_Load = "False"
Incremental_Displacement = "False"

#PostProcess Data
#####################################

nodal_results=["DISPLACEMENT","REACTION","PRESSURE"]
gauss_points_results=["GREEN_LAGRANGE_STRAIN_TENSOR","CAUCHY_STRESS_TENSOR","FORCE"]

# GiD output configuration
class GidOutputConfiguration:
    GiDPostMode = "Binary"
    GiDWriteMeshFlag = True
    GiDWriteConditionsFlag = True
    GiDWriteParticlesFlag = False
    GiDMultiFileFlag = "Multiples"

GiDWriteFrequency = 1
WriteResults = "PreMeshing"
echo_level = 1

# graph_options
PlotGraphs = "False"
PlotFrequency = 10

# list options
PrintLists = "False"
number_of_lists = 1
file_list = []
file_list.append(10)

# restart options
SaveRestart = "False"
RestartFrequency = 0
LoadRestart = "False"
Restart_Step = 0

# Declare Python Variables

problem_name = 'cooks_membrane_test_2D'
problem_path = '/home/jmaria/ModelsGID/PfemMechanicsApplication/cooks_membrane_test_2D.gid'
kratos_path = '/home/jmaria'

