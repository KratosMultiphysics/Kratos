from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
domain_size = 3

Convergence_Tolerance = 1.0
Absolute_Tolerance = 1.0
Max_Iter = 1
Convergence_Criteria = "Displacement_Criteria"

Residual_Convergence_Tolerance = 1.0E-3
Residual_Absolute_Tolerance = 1.0E-6
Displacement_Convergence_Tolerance = 1.0E-6
Displacement_Absolute_Tolerance = 1.0E-9
Solution_method = "Newton-Raphson"
SolverType = "StaticSolver"
LinearSolver = "SuperLUSolver"
FindNodalNeighbours = "False"
FindElementalNeighbours = "False"
Rotational_Dofs = "True"

nodal_results = ["DISPLACEMENT", "REACTION"]
gauss_points_results = ["GREEN_LAGRANGE_STRAIN_TENSOR", "ROTATION", "PK2_STRESS_TENSOR", "MOMENT"]
GiDPostMode = "Binary"
GiDWriteMeshFlag = True
GiDWriteConditionsFlag = True
GiDMultiFileFlag = "Single"

problem_name = "scordelisStructuralAnalysis"
problem_path = "/home/rrossi/esempi/scordelis.gid"
