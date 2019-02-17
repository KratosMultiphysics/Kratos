from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
domain_size = 2

Dt = 0.00001
max_time = 1
nsteps = 400
FindNodalNeighbours = "False"
FindElementalNeighbours = "True"
Rotational_Dofs = "False"
SolverType = "StaticSolver"
LinearSolver = "SkylineLUFactorization"
Linear_Solver_Tolerance = 1E-9
Linear_Solver_Max_Iteration = 5000
Convergence_Criteria = "Displacement_Criteria"
Convergence_Tolerance = 1E-6
Absolute_Tolerance = 1E-9
Max_Iter = 30
# Declare Python Variables

problem_name = "plasticJ2"
problem_path = "."
kratos_path = "../../../.."
