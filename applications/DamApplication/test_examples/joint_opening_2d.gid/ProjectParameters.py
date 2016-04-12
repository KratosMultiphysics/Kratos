## General Data ---------------------------------------------------------------------------------------------------------------

domain_size = 2
NumberofThreads = 1
time_scale = "Seconds"
evolution_type = "Exact" 
delta_time = 0.01
ending_time = 0.02


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
    residual_relative_tolerance = 1.0E-4
    max_iteration = 15
    linear_solver = "Direct"
    direct_solver = "Super_LU"
    iterative_solver = "AMGCL"
    compute_reactions = True


## Boundary Conditions Data ---------------------------------------------------------------------------------------------------

class ConditionsOptions:
    Imposed_Displacement = "Unmodified"
    Imposed_Temperature ="Unmodified"
    Imposed_PointLoad = "Unmodified"
    Imposed_LineLoad = "Unmodified"
    Imposed_SurfaceLoad = "Unmodified"
    Imposed_NormalLoad = "Unmodified"
    Imposed_WaterLoad = "Unmodified"
    Imposed_Bofang_Temperature = "Unmodified"

    
## PostProcess Data -----------------------------------------------------------------------------------------------------------

nodal_results = ["DISPLACEMENT","REACTION","TEMPERATURE","NO_RESULT","LINE_LOAD","NO_RESULT","NO_RESULT"]
gauss_points_results=["NO_RESULT","NO_RESULT","CAUCHY_STRESS_TENSOR","MECHANICAL_STRESS_TENSOR","THERMAL_STRESS_TENSOR","NO_RESULT","DAMAGE_VARIABLE"]

class GidOutputConfiguration:
    GiDWriteMeshFlag = True
    GiDWriteConditionsFlag = False
    GiDPostFiles = "Single"
    GiDPostMode = "Binary"
    GiDWriteFrequency = 0.01


## Problem Data ---------------------------------------------------------------------------------------------------------------

problem_name = 'joint_opening_2d'
problem_path = '/path/to/kratos/applications/DamApplication/test_examples/joint_opening_2d.gid'

