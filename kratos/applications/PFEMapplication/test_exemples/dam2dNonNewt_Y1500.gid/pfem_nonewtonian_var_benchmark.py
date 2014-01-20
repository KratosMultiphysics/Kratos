from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
domain_size = 2

Density = 1000.0
Viscosity = 0.000001
Gravity_X = 0.0
Gravity_Y = -9.806
Gravity_Z = 0.0
YieldStress = 1500.0
max_x = 2.5
max_y = 3.0
max_z = 1.0
min_x = -1.0
min_y = -1.0
min_z = -1.0
use_oss = 0
dynamic_tau = 1
nsteps = 150
output_Dt = 0.5
min_dt = 0.01
max_dt = 0.1
safety_factor = 0.1
SolverType = "monolithic_solver_lagrangian"
laplacian_form = 3
# Declare Python Variables

problem_name = "dam2dNonNewt_Y1500"
problem_path = "/home/antonia/kratos/applications/PFEMapplication/test_exemples/dam2dNonNewt_Y1500.gid"
kratos_path = "/home/antonia/kratos"
