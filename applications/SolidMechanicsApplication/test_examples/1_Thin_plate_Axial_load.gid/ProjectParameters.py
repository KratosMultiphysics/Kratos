domain_size = 2

#Problem Data
#################################################

ProblemType = "Mechanical"
NumberofThreads = 2
Solution_method = "Newton-Raphson"
SolverType = "StaticSolver"
time_step = 1
end_time    = 1

#Solver Data
#################################################

class SolverSettings:
    solver_type = "mechanical_solver"
    domain_size = 2
    echo_level  = 0

    max_delta_time  = time_step
    time_integration_method = "Implicit"
    explicit_integration_scheme = "CentralDifferences"
    time_step_prediction_level  = "Automatic"
    rayleigh_damping = False

    RotationDofs = False
    PressureDofs = False
    ReformDofSetAtEachStep = False
    LineSearch = False
    Implex = False
    ComputeReactions = True
    ComputeContactForces = False
    scheme_type = "StaticSolver"
    convergence_criterion = "Residual_criteria"
    displacement_relative_tolerance = 1.0E-4
    displacement_absolute_tolerance = 1.0E-9
    residual_relative_tolerance = 1.0E-4
    residual_absolute_tolerance = 1.0E-9
    max_iteration = 10
    class linear_solver_config:
        solver_type = "Super LU"
        scaling = False


#Constraints Data
#################################################

Incremental_Load = "False"
Incremental_Displacement = "False"

#PostProcess Data
#################################################

nodal_results=["DISPLACEMENT","REACTION"]
gauss_points_results=[]

# GiD output configuration
class GidOutputConfiguration:
    GiDPostMode = "Binary"
    GiDWriteMeshFlag = False
    GiDWriteConditionsFlag = False
    GiDWriteParticlesFlag = False
    GiDMultiFileFlag = "Single"

GiDWriteFrequency = 0.01
WriteResults = "PreMeshing"
echo_level = 0

# graph_options
PlotGraphs = "False"
PlotFrequency = 0 

# list options
PrintLists = "True"
file_list = [] 

# restart options
SaveRestart = False
RestartFrequency = 0
LoadRestart = False
Restart_Step = 0

problem_name="1_Thin_plate_Axial_loadStructuralAnalysis"
problem_path="/home/flo/Escritorio/Josep Maria/Gid/1_Thin_plate_Axial_load.gid"
