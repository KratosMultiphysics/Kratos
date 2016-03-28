## General Data ---------------------------------------------------------------------------------------------------------------
NumberofThreads = 4
Start_time = 0
Dt = 1.0
max_time = 1.0
output_time = 1.0
nsteps =  1
ProblemType = "Mechanical"

## Solver Data ----------------------------------------------------------------------------------------------------------------

linear_solver = "Direct"
direct_solver = "Super_LU"
iterative_solver = "AMGCL"

class SolidSolverConfiguration:
    solver_type = "mechanical_solver"
    echo_level  = 1
    domain_size = 3 
    
    time_integration_method = "Implicit"
    explicit_integration_scheme = "CentralDifferences"
    time_step_prediction_level  = 0

    LineSearch =  False
    Implex =  False
    ComputeReactions = False 
    ComputeContactForces = False 
    scheme_type = "StaticSolver"
    convergence_criterion = "Residual_criteria" 
    displacement_relative_tolerance = 1.0E-4
    displacement_absolute_tolerance =1.0E-4
    residual_relative_tolerance =  1.0E-4
    residual_absolute_tolerance =  1.0E-4
    max_iteration = 30
    if linear_solver == "Direct":
        class linear_solver_config:
            solver_type = direct_solver
            scaling = False
    else:
        class linear_solver_config:
            solver_type = iterative_solver
            max_iteration = 100
            tolerance = 1.0E-4

## Boundary Conditions Data -----------------------------------------------------------------------------------------------------------

class ConditionsOptions:
    Imposed_Displacement = "Constant"
    Imposed_PointLoad = "Constant"
    
## Problem Data ---------------------------------------------------------------------------------------------------------------
kratos_path="~/bin/kratos"
problem_name = 'sphere_test'
problem_path = '/home/vicente/Descargas/sphere_test.gid'

