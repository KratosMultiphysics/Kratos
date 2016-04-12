## General Data ---------------------------------------------------------------------------------------------------------------

domain_size = 3
NumberofThreads = 1
time_scale = "Hours"
evolution_type = "Exact" 
delta_time = 1
ending_time = 5


## Solver Data ----------------------------------------------------------------------------------------------------------------

class DiffusionSolverSettings:
    unknown_variable = "TEMPERATURE"
    diffusion_variable = "CONDUCTIVITY"
    specific_heat_variable = "SPECIFIC_HEAT"
    density_variable = "DENSITY"
    
class MechanicalSolverSettings:
    analysis_type = "Quasi-Static"
    strategy_type = "Newton-Raphson"
    dofs_relative_tolerance = 1.0E-4
    residual_relative_tolerance = 1.0E-3
    max_iteration = 10
    linear_solver = "Iterative"
    direct_solver = "Super_LU"
    iterative_solver = "AMGCL"
    compute_reactions = False


## Boundary Conditions Data ---------------------------------------------------------------------------------------------------

class ConditionsOptions:
    Imposed_Displacement = "Unmodified"
    Imposed_Temperature ="Unmodified"
    Imposed_PointLoad = "Unmodified"
    Imposed_LineLoad = "Unmodified"
    Imposed_SurfaceLoad = "Unmodified"
    Imposed_NormalLoad = "Unmodified"
    Imposed_WaterLoad = "Table_Evolution_Data"
    Imposed_Bofang_Temperature = "Unmodified"

    
## PostProcess Data -----------------------------------------------------------------------------------------------------------

nodal_results = ["DISPLACEMENT","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT","NEGATIVE_FACE_PRESSURE"]
gauss_points_results=["GREEN_LAGRANGE_STRAIN_TENSOR","NO_RESULT","CAUCHY_STRESS_TENSOR","NO_RESULT","NO_RESULT","VON_MISES_STRESS","NO_RESULT"]

class GidOutputConfiguration:
    GiDWriteMeshFlag = True
    GiDWriteConditionsFlag = False
    GiDPostFiles = "Single"
    GiDPostMode = "Binary"
    GiDWriteFrequency = 1


## Problem Data ---------------------------------------------------------------------------------------------------------------

problem_name = 'upliftpressure_3D'
problem_path = '/path/to/kratos/applications/DamApplication/test_examples/upliftpressure_3D.gid'

