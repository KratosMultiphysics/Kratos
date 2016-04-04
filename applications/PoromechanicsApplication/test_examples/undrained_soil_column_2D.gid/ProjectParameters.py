## General Data ---------------------------------------------------------------------------------------------------------------

domain_size = 2
NumberofThreads = 1
delta_time = 0.02
ending_time = 2
fic_stabilization = False


## Solver Data ----------------------------------------------------------------------------------------------------------------

class SolverSettings:
    analysis_type = "Quasi-Static"
    strategy_type = "Newton-Raphson"
    dofs_relative_tolerance = 1.0E-4
    residual_relative_tolerance = 1.0E-4
    max_iteration = 15
    linear_solver = "Direct"
    direct_solver = "Super_LU"
    iterative_solver = "AMGCL"
    compute_reactions = False


## Boundary Conditions Data ---------------------------------------------------------------------------------------------------

class ConditionsOptions:
    Imposed_Displacement = "Unmodified"
    Imposed_Pressure = "Unmodified"
    Imposed_Force = "Unmodified"
    Imposed_FaceLoad = "Table_Interpolation"
    Imposed_NormalLoad = "Unmodified"
    Imposed_TangentialLoad = "Unmodified"
    Imposed_NormalFluidFlux = "Unmodified"


## PostProcess Data -----------------------------------------------------------------------------------------------------------

nodal_results = ["DISPLACEMENT","WATER_PRESSURE","NO_RESULT","NO_RESULT","NO_RESULT","FACE_LOAD","NO_RESULT","NO_RESULT","NO_RESULT"]
gauss_points_results=["NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT","NO_RESULT"]

class GidOutputConfiguration:
    GiDWriteMeshFlag = True
    GiDWriteConditionsFlag = False
    GiDPostFiles = "Single"
    GiDPostMode = "Binary"
    GiDWriteFrequency = 0.02


## Problem Data ---------------------------------------------------------------------------------------------------------------

problem_name = 'undrained_soil_column_2D'
problem_path = '/path/to/kratos/applications/PoromechanicsApplication/test_examples/undrained_soil_column_2D.gid'

