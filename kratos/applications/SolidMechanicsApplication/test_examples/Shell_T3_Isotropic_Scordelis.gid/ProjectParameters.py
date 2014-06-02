domain_size = 3

#Problem Data
#################################################

ProblemType = "Mechanical"
NumberofThreads = 2
Solution_method = "Newton-Raphson"
SolverType = "StaticSolver"
time_step = 1
nsteps    = 1

#Solver Data
#################################################

class SolverSettings:
    solver_type = "mechanical_solver"
    domain_size = 3
    echo_level  = 0
    RotationDofs = True
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
gauss_points_results=["MOMENT","FORCE"]

# GiD output configuration
class GidOutputConfiguration:
    GiDPostMode = "Binary"
    GiDWriteMeshFlag = True
    GiDWriteConditionsFlag = True
    GiDWriteParticlesFlag = False
    GiDMultiFileFlag = "Single"

GiDWriteFrequency = 0.01
WriteResults = "PreMeshing"
echo_level = 1

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

problem_name="Shell_T3_Isotropic_ScordelisStructuralAnalysis"
problem_path="/home/jmaria/kratos/applications/SolidMechanicsApplication/test_examples/Shell_T3_Isotropic_Scordelis.gid"
