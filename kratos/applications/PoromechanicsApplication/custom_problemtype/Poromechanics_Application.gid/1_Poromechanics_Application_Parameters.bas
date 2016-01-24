## General Data ---------------------------------------------------------------------------------------------------------------

domain_size = *GenData(Domain_Size,INT)
NumberofThreads = *GenData(Number_of_threads,INT)
delta_time = *GenData(Delta_Time)
ending_time = *GenData(Ending_Time)
fic_stabilization = *GenData(FIC_Stabilization)


## Solver Data ----------------------------------------------------------------------------------------------------------------

class SolverSettings:
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
    Imposed_Pressure = "*GenData(Imposed_Pressure)"
    Imposed_PointLoad = "*GenData(Imposed_PointLoad)"
    Imposed_LineLoad = "*GenData(Imposed_LineLoad)"
    Imposed_SurfaceLoad = "*GenData(Imposed_SurfaceLoad)"
    Imposed_NormalLoad = "*GenData(Imposed_NormalLoad)"
    Imposed_TangentialLoad = "*GenData(Imposed_TangentialLoad)"
    Imposed_NormalFluidFlux = "*GenData(Imposed_NormalFluidFlux)"


## PostProcess Data -----------------------------------------------------------------------------------------------------------

nodal_results = ["*GenData(Nodal_results_1)","*GenData(Nodal_results_2)","*GenData(Nodal_results_3)","*GenData(Nodal_results_4)","*GenData(Nodal_results_5)","*GenData(Nodal_results_6)","*GenData(Nodal_results_7)","*GenData(Nodal_results_8)","*GenData(Nodal_results_9)","*GenData(Nodal_results_10)"]
gauss_points_results=["*GenData(Gauss_points_results_1)","*GenData(Gauss_points_results_2)","*GenData(Gauss_points_results_3)","*GenData(Gauss_points_results_4)"]

class GidOutputConfiguration:
    GiDWriteMeshFlag = *GenData(Write_deformed_mesh)
    GiDWriteConditionsFlag = *GenData(Write_conditions)
    GiDPostFiles = "*GenData(GiD_post_files)"
    GiDPostMode = "*GenData(GiD_post_mode)"
    GiDWriteFrequency = *GenData(Write_Frequency)


## Problem Data ---------------------------------------------------------------------------------------------------------------

