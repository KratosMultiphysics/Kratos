domain_size = 3

#Problem Data
#################################################

ProblemType = "Mechanical"
NumberofThreads = 8
Solution_method = "Newton-Raphson"
SolverType = "QuasiStaticSolver"
time_step = 0.05
nsteps    = 20

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
    scheme_type = "QuasiStaticSolver"
    convergence_criterion = "Residual_criteria"
    displacement_relative_tolerance = 1.0E-6
    displacement_absolute_tolerance = 1.0E-9
    residual_relative_tolerance = 1.0E-6
    residual_absolute_tolerance = 1.0E-6
    max_iteration = 30
    class linear_solver_config:
        solver_type = "Super LU"
        scaling = False


#Constraints Data
#################################################

Incremental_Load = "True"
Incremental_Displacement = "True"

#PostProcess Data
#################################################

nodal_results=["DISPLACEMENT","REACTION"]
gauss_points_results=["SHELL_FORCE","SHELL_FORCE_GLOBAL","SHELL_MOMENT","SHELL_STRAIN","SHELL_STRAIN_GLOBAL","MATERIAL_ORIENTATION_DX","MATERIAL_ORIENTATION_DY","MATERIAL_ORIENTATION_DZ"]

# GiD output configuration
class GidOutputConfiguration:
    GiDPostMode = "Binary"
    GiDWriteMeshFlag = False
    GiDWriteConditionsFlag = False
    GiDWriteParticlesFlag = False
    GiDMultiFileFlag = "Single"

GiDWriteFrequency = 0.05
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

problem_name="Shell_T3_Thin__BendingRollUpStructuralAnalysis"
problem_path="/home/kratos/applications/SolidMechanicsApplication/test_examples/Shell_T3_Thin__BendingRollUp.gid"
