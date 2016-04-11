## General Data ---------------------------------------------------------------------------------------------------------------

domain_size = *GenData(Domain_Size,INT)
NumberofThreads = *GenData(Number_of_threads,INT)
time_scale = "*GenData(Time_Scale)"
evolution_type = "*GenData(Evolution_Type)" 
delta_time = *GenData(Delta_Time)
ending_time = *GenData(Ending_Time)


## Solver Data ----------------------------------------------------------------------------------------------------------------

class DiffusionSolverSettings:
    unknown_variable = "TEMPERATURE"
    diffusion_variable = "CONDUCTIVITY"
    specific_heat_variable = "SPECIFIC_HEAT"
    density_variable = "DENSITY"
    
class MechanicalSolverSettings:
    analysis_type = "*GenData(Analysis_Type)"
    strategy_type = "*GenData(Strategy_Type)"
    dofs_relative_tolerance = *GenData(DOFs_Relative_Tolerance)
    residual_relative_tolerance = *GenData(Residual_Relative_Tolerance)
    max_iteration = *GenData(Max_Iterations,INT)
    linear_solver = "*GenData(Linear_Solver)"
    direct_solver = "*GenData(Direct_Solver_Type)"
    iterative_solver = "*GenData(Iterative_Solver_Type)"
    compute_reactions = *GenData(Write_Reactions)


## Boundary Conditions Data ---------------------------------------------------------------------------------------------------

class ConditionsOptions:
    Imposed_Displacement = "*GenData(Imposed_Displacement)"
    Imposed_Temperature ="*GenData(Imposed_Temperature)"
    Imposed_PointLoad = "*GenData(Imposed_PointLoad)"
    Imposed_LineLoad = "*GenData(Imposed_LineLoad)"
    Imposed_SurfaceLoad = "*GenData(Imposed_SurfaceLoad)"
    Imposed_NormalLoad = "*GenData(Imposed_NormalLoad)"
    Imposed_WaterLoad = "*GenData(Imposed_WaterLoad)"
    Imposed_Bofang_Temperature = "*GenData(Imposed_Bofang_Temperature)"

    
## PostProcess Data -----------------------------------------------------------------------------------------------------------

nodal_results = ["*GenData(Nodal_results_1)","*GenData(Nodal_results_2)","*GenData(Nodal_results_3)","*GenData(Nodal_results_4)","*GenData(Nodal_results_5)","*GenData(Nodal_results_6)","*GenData(Nodal_results_7)"]
gauss_points_results=["*GenData(Gauss_points_results_1)","*GenData(Gauss_points_results_2)","*GenData(Gauss_points_results_3)","*GenData(Gauss_points_results_4)","*GenData(Gauss_points_results_5)","*GenData(Gauss_points_results_6)","*GenData(Gauss_points_results_7)"]

class GidOutputConfiguration:
    GiDWriteMeshFlag = *GenData(Write_deformed_mesh)
    GiDWriteConditionsFlag = *GenData(Write_conditions)
    GiDPostFiles = "*GenData(GiD_post_files)"
    GiDPostMode = "*GenData(GiD_post_mode)"
    GiDWriteFrequency = *GenData(Write_Frequency)


## Problem Data ---------------------------------------------------------------------------------------------------------------

