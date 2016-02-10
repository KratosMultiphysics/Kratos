## General Data ---------------------------------------------------------------------------------------------------------------

domain_size = 3
plane_state = "Plane_Stress"
NumberofThreads = 1
time_scale = "Seconds"
evolution_type = "Exact" 
delta_time = 0.002
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
    residual_relative_tolerance = 1.0E-3
    max_iteration = 15
    linear_solver = "Iterative"
    direct_solver = "Super_LU"
    iterative_solver = "AMGCL"
    compute_reactions = True


## Boundary Conditions Data ---------------------------------------------------------------------------------------------------

class ConditionsOptions:
    Imposed_Displacement = "Constant"
    Imposed_PointLoad = "Constant"
    Imposed_LineLoad = "Constant"
    Imposed_SurfaceLoad = "Constant"
    Imposed_NormalLoad = "Constant"
    Imposed_TangentialLoad = "Constant"
    Imposed_Temperature ="Constant"
    
## PostProcess Data -----------------------------------------------------------------------------------------------------------

nodal_results = ["DISPLACEMENT","REACTION","TEMPERATURE","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT"]
gauss_points_results=["NO_RESULT","NO_RESULT","CAUCHY_STRESS_TENSOR","MECHANICAL_STRESS_TENSOR","THERMAL_STRESS_TENSOR","VON_MISES_STRESS"]

class GidOutputConfiguration:
    GiDWriteMeshFlag = True
    GiDWriteConditionsFlag = False
    GiDPostFiles = "Single"
    GiDPostMode = "Binary"
    GiDWriteFrequency = 0.002


## Problem Data ---------------------------------------------------------------------------------------------------------------

problem_name = '/path/to/kratos/applications/DamApplication/test_examples/diffusion_beam_3d.gid/diffusion_beam_3d'
problem_path = '/path/to/kratos/applications/DamApplication/test_examples/diffusion_beam_3d.gid'

