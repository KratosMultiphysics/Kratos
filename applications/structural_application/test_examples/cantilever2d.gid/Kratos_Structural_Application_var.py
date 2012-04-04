domain_size = 2

Dt = 0.01
max_time = 1
nsteps = 5
FindNodalNeighbours = "False"
FindElementalNeighbours = "False"
Rotational_Dofs = "False"
SolverType = "StaticSolver"
TimeIntegration = "Bossak_Newmark"
LinearSolver = "SkylineLUFactorization"
Linear_Solver_Tolerance = 1E-9
Linear_Solver_Max_Iteration = 5000
Convergence_Criteria = "Displacement_Criteria"
Convergence_Tolerance = 1E-6
Absolute_Tolerance = 1E-9
Max_Iter = 30
# Declare Python Variables

problem_name="cantilever2d"
