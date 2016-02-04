## General Data ---------------------------------------------------------------------------------------------------------------

domain_size = 2
plane_state = "Plane_Strain"
NumberofThreads = 1
time_scale = "Hours"
evolution_type = "Interpolation" 
delta_time = 0.5
ending_time = 3


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
    linear_solver = "Direct"
    direct_solver = "Super_LU"
    iterative_solver = "AMGCL"
    compute_reactions = True


## Boundary Conditions Data ---------------------------------------------------------------------------------------------------

class ConditionsOptions:
    Imposed_Displacement = "Constant"
    Imposed_PointLoad = "Constant"
    Imposed_LineLoad = "Constant"
    Imposed_SurfaceLoad = "Constant"
    Imposed_NormalLoad = "Table"
    Imposed_TangentialLoad = "Constant"
    Imposed_Temperature ="Table"
    
## PostProcess Data -----------------------------------------------------------------------------------------------------------

nodal_results = ["DISPLACEMENT","REACTION","TEMPERATURE","NO_RESULT","NO_RESULT","NO_RESULT","NORMAL_CONTACT_STRESS","NO_RESULT"]
gauss_points_results=["GREEN_LAGRANGE_STRAIN_TENSOR","THERMAL_STRAIN_TENSOR","CAUCHY_STRESS_TENSOR","MECHANICAL_STRESS_TENSOR","THERMAL_STRESS_TENSOR","VON_MISES_STRESS"]

class GidOutputConfiguration:
    GiDWriteMeshFlag = True
    GiDWriteConditionsFlag = False
    GiDPostFiles = "Single"
    GiDPostMode = "Binary"
    GiDWriteFrequency = 0.5


## Problem Data ---------------------------------------------------------------------------------------------------------------

problem_name = 'thermomechanical_dam_2d'
problem_path = '/path/to/kratos/applications/DamApplication/test_examples/thermomechanical_dam_2d.gid'

