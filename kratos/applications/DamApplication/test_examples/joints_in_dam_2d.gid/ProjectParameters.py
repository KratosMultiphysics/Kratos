## General Data ---------------------------------------------------------------------------------------------------------------

domain_size = 2
NumberofThreads = 1
time_scale = "Seconds"
evolution_type = "Exact" 
delta_time = 1
ending_time = 10


## Solver Data ----------------------------------------------------------------------------------------------------------------

class DiffusionSolverSettings:
    unknown_variable = "TEMPERATURE"
    diffusion_variable = "CONDUCTIVITY"
    specific_heat_variable = "SPECIFIC_HEAT"
    density_variable = "DENSITY"
    
class MechanicalSolverSettings:
    analysis_type = "Quasi-Static"
    strategy_type = "Newton-Raphson"
    convergence_criterion = "Displacement_criterion"
    displacement_rel_tol = 1.0E-4
    displacement_abs_tol = 1.0E-8
    residual_rel_tol = 1.0E-4
    residual_abs_tol = 1.0E-8
    echo_level = 1
    max_iteration = 30
    linear_solver = "Direct"
    direct_solver = "Super_LU"
    iterative_solver = "BICGSTAB"
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

nodal_results = ["DISPLACEMENT","NO_RESULT","TEMPERATURE","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT"]
gauss_points_results=["NO_RESULT","NO_RESULT","CAUCHY_STRESS_TENSOR","MECHANICAL_STRESS_TENSOR","THERMAL_STRESS_TENSOR","NO_RESULT","NO_RESULT"]

class GidOutputConfiguration:
    GiDWriteMeshFlag = True
    GiDWriteConditionsFlag = False
    GiDPostFiles = "Single"
    GiDPostMode = "Binary"
    GiDWriteFrequency = 1


## Problem Data ---------------------------------------------------------------------------------------------------------------

problem_name = 'joints_in_dam_2d'
problem_path = '/path/to/kratos/applications/DamApplication/test_examples/joints_in_dam_2d.gid'

