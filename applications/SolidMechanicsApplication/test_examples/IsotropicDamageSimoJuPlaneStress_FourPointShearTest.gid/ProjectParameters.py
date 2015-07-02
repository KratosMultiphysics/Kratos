domain_size = 2

#Problem Data
#################################################

ProblemType = "Mechanical"
NumberofThreads = 2
Solution_method = "Newton-Raphson"
SolverType = "QuasiStaticSolver"
time_step = 0.01
end_time    = 2.4

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

    RotationDofs = False
    PressureDofs = False
    ReformDofSetAtEachStep = False
    LineSearch = False
    Implex = False
    ComputeReactions = True
    ComputeContactForces = False
    scheme_type = "QuasiStaticSolver"
    convergence_criterion = "And_criteria"
    displacement_relative_tolerance = 1.0E-6
    displacement_absolute_tolerance = 1.0E-9
    residual_relative_tolerance = 1.0E-5
    residual_absolute_tolerance = 1.0E-7
    max_iteration = 30
    class linear_solver_config:
        solver_type = "Super LU"
        scaling = False


#Constraints Data
#################################################

Incremental_Load = "True"
Incremental_Displacement = "False"

#PostProcess Data
#################################################

nodal_results=["DISPLACEMENT","LINE_LOAD"]
gauss_points_results=["GREEN_LAGRANGE_STRAIN_TENSOR","CAUCHY_STRESS_TENSOR","DAMAGE_VARIABLE"]

# GiD output configuration
class GidOutputConfiguration:
    GiDPostMode = "Binary"
    GiDWriteMeshFlag = False
    GiDWriteConditionsFlag = False
    GiDWriteParticlesFlag = False
    GiDMultiFileFlag = "Single"

GiDWriteFrequency = 0.02
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

problem_name="IsotropicDamageSimoJuPlaneStress_FourPointShearTestStructuralAnalysis"
problem_path="/home/ipouplana/kratos/applications/SolidMechanicsApplication/test_examples/IsotropicDamageSimoJuPlaneStress_FourPointShearTest.gid"
